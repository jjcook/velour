//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// pregraph_distribution.cpp
//

// TODO: assert that node and sequence(read) count for pregraph file is less than LONG_MAX for Velvet


// TODO:
//   # special case, 1 partition, to still do read filtering
//

#include "types.h"

static unsigned *vid_to_partition_map = NULL;

static const unsigned MAXPARTS = 64;

static FILE *pregraphFile[MAXPARTS] = {NULL};
static FILE *readsFile[MAXPARTS] = {NULL};

static uint64_t nodeCount[MAXPARTS] = {0};
static uint64_t readCount[MAXPARTS] = {0};

static unsigned current_partition;
static uint64_t hashmask;
static uint64_t *hashtable;

static inline void hash_insert(uint64_t bucket, unsigned value)
{
    assert( value < 64 );
    assert( bucket <= hashmask );
    hashtable[bucket] |= (1 << value);
}

static inline uint64_t hash_fetch(uint64_t bucket)
{
    assert( bucket <= hashmask );
    return hashtable[bucket];
}

static inline uint64_t hash(const Kmer kmer) {
    return ((kmer ^ (kmer >> 17) ^ (kmer >> 37) ^ (kmer << 11)) & hashmask);
}

namespace {
struct lambda_fillhash {
    lambda_fillhash() {}
    void operator()(Kmer k, Kmer canon_k, bool new_sense_reversed, Nucleotide base, Nucleotide last_base, unsigned kmer_length, bool isFirstKmer, bool isLastKmer)
    {
        hash_insert(hash(canon_k), current_partition);
    }
};

struct lambda_hashread {
    lambda_hashread() : readmask(0) {}
    void operator()(Kmer k, Kmer canon_k, bool new_sense_reversed, Nucleotide base, Nucleotide last_base, unsigned kmer_length, bool isFirstKmer, bool isLastKmer)
    {
        readmask |= hash_fetch(hash(canon_k));
    }
    uint64_t readmask;
};
} // namespace: anonymous


static void emitPartitionSequence(unsigned partition, unsigned category, char *header, Sequence *sequence)
{
    FILE *file = readsFile[partition];
    fprintf(file, "%s\t%"PRIu64"\t%u\n", header, ++(readCount[partition]), category);
    sequence->PrintToFile(file);
    fputc('\n', file);
}
 
static void distribute_reads(unsigned partition_count, unsigned kmer_length, FILE *allreads)
{
    Sequence_StackAllocated memory1;
    Sequence_StackAllocated memory2;

    unsigned category1;
    char header1[MAX_READ_HEADER_LENGTH];
    Sequence *sequence1 = new (&memory1) Sequence(Sequence::MAX_BASES);

    unsigned category2;
    char header2[MAX_READ_HEADER_LENGTH];
    Sequence *sequence2 = new (&memory2) Sequence(Sequence::MAX_BASES);

    char ascii_sequence[MAX_READ_LENGTH];
    unsigned ascii_offset = 0;

    char line[MAX_READ_LENGTH];
    char *retval = fgets(line, MAX_READ_LENGTH, allreads);
    while (retval != NULL) {
        assert( line[0] == '>' );
        sscanf(line, "%[^\t]\t%*u\t%u", header1, &category1);
        ascii_offset = 0;
        ascii_sequence[0] = '\0';
        while ((retval = fgets(line, MAX_READ_LENGTH, allreads)) != NULL) {
            if (line[0] == '>') break;
            sscanf(line, "%s", &ascii_sequence[ascii_offset]);
            ascii_offset += strlen(&ascii_sequence[ascii_offset]);
            assert( ascii_offset < MAX_READ_LENGTH );
        }
        ascii_sequence[ascii_offset] = '\0';
        sequence1->InitializeWithString(ascii_sequence);

        // get paired read
        bool ispaired = (category1 % 2 == 1) && (category1 < 2*CATEGORIES);
        if (ispaired) {
            assert( line[0] == '>' );
            sscanf(line, "%[^\t]\t%*u\t%u", header2, &category2);
            assert( category1 == category2 );
            ascii_offset = 0;
            ascii_sequence[0] = '\0';
            while ((retval = fgets(line, MAX_READ_LENGTH, allreads)) != NULL) {
                if (line[0] == '>') break;
                sscanf(line, "%s", &ascii_sequence[ascii_offset]);
                ascii_offset += strlen(&ascii_sequence[ascii_offset]);
                assert( ascii_offset < MAX_READ_LENGTH );
            }
            ascii_sequence[ascii_offset] = '\0';
            sequence2->InitializeWithString(ascii_sequence);
        }

        // compute partition mask for the read(s)
        lambda_hashread functor_readmask1;
        sequence_process_kmers(sequence1, kmer_length, functor_readmask1);

        lambda_hashread functor_readmask2;
        if (ispaired) {
            sequence_process_kmers(sequence2, kmer_length, functor_readmask2);
        }

        uint64_t pairmask = functor_readmask1.readmask & functor_readmask2.readmask; // partitions both reads map to

        // distribute read(s) to each partition
        for (unsigned i=0 ; i < partition_count ; ++i) {
            if (ispaired && ((pairmask >> i) & 0x1)) { // pair matches
                assert( category1 == category2 );
                emitPartitionSequence(i, category1, header1, sequence1);
                emitPartitionSequence(i, category2, header2, sequence2);
            } else {
                if ((functor_readmask1.readmask >> i) & 0x1) {
                    unsigned nonpaired_category = category1 & (~0x1);
                    assert( nonpaired_category < 100 );
                    emitPartitionSequence(i, nonpaired_category, header1, sequence1);
                }
                if ((functor_readmask2.readmask >> i) & 0x1) {
                    unsigned nonpaired_category = category2 & (~0x1);
                    emitPartitionSequence(i, nonpaired_category, header2, sequence2);
                }
            }
        }
    }
    assert( feof(allreads) );
}


void runPregraphDistribution(file_object_vector *file_objects, SeqGraph *sgraph)
{
    if (g__PGDIST_PARTITIONS > MAXPARTS) {
        fprintf(stderr, "ERROR: PGDIST partitions must be less than %d if read filtering.  Exiting...\n", MAXPARTS);
        exit(EXIT_FAILURE);
    }

    // initialize the read count for pregraph emittal
    {
        char filename[PATH_MAX+1];
        sprintf(filename, "%s/common.sequences", g__WORK_BASE_DIRECTORY);
        FILE *f = fopen(filename, "r");
        if (f == NULL) {
            fprintf(stderr, "ERROR: fopen()\n");
            assert(false);
            exit(EXIT_FAILURE);
        }
        fscanf(f, "%"SCNu64"", &g__READ_COUNT);
        fclose(f);
    }

    // load the pregraph quilt, with hypergraph vids
    for(file_object_vector::iterator itr = file_objects->begin() ; itr != file_objects->end() ; ++itr) {
        switch (itr->filetype) {
        case QUILT:
            sg_load_quilt(sgraph, itr->filename);
            break;
        default:
            fprintf(stderr, "ERROR: Cannot use %s file %s as pregraph distribution input.  Exiting...\n",
                FILE_TYPES[itr->filetype], itr->filename);
            exit(EXIT_FAILURE);
        }
    }

    // open the metis file to determine vid count, and allocate map
    uint64_t vid_count;
    {
        char filename[PATH_MAX+1];
        sprintf(filename, "%s/PreGraph.metis", g__WORK_BASE_DIRECTORY);
        FILE *f = fopen(filename, "r");
        if (f == NULL) {
            fprintf(stderr, "ERROR: fopen()\n");
            assert(false);
            exit(EXIT_FAILURE);
        }
        fscanf(f, "%"SCNu64"", &vid_count);
        fclose(f);
    }
    vid_to_partition_map = (unsigned *) calloc(vid_count, sizeof(unsigned));
    assert( vid_to_partition_map != NULL );

    // load the vid --> part mapping
    {
        char filename[PATH_MAX+1];
        sprintf(filename, "%s/PreGraph.metis.part.%u", g__WORK_BASE_DIRECTORY, g__PGDIST_PARTITIONS);
        FILE *f = fopen(filename, "r");
        if (f == NULL) {
            fprintf(stderr, "ERROR: fopen()\n");
            assert(false);
            exit(EXIT_FAILURE);
        }

        for (uint64_t i=1; i <= vid_count; ++i) {
            fscanf(f, "%u\n", &vid_to_partition_map[i]);
        }
        fclose(f);
    }

    // apply vid --> part map to nodes via iterator
    struct lambda1 {
        static void apply_map(SeqNode *node)
        {
            node->right_vid = vid_to_partition_map[abs(node->right_vid)];
            node->left_vid = vid_to_partition_map[abs(node->left_vid)];
        }
    };
    sg_for_each(sgraph, lambda1::apply_map);

    // deallocate vid map
    free(vid_to_partition_map);

    //
    // initialize the pregraph files
    //
    for (unsigned i=0 ; i < g__PGDIST_PARTITIONS; ++i) {
        char pregraph_filename[PATH_MAX+1];
        sprintf(pregraph_filename, "%s/%u-PreGraph", g__WORK_BASE_DIRECTORY, i);
        int pregraph_filedes = open(pregraph_filename, O_CREAT | O_TRUNC | O_WRONLY, 0666);
        pregraphFile[i] = fdopen(pregraph_filedes, "w");

        nodeCount[i] = 0;
        readCount[i] = 0;

        // PREGRAPH File Header: node count, read count, kmer length, 1
        //   note: using hex so can fixup node and sequence count later
        fprintf(pregraphFile[i], "0x%08lx\t0x%08lx\t%u\t%hd\n", 0UL, (long) g__READ_COUNT, g__FULLKMER_LENGTH, 1);
    }

    //  - distribute pregraph nodes via iterator
    struct lambda2 {
        static void emit_pregraph_sequence(SeqNode *node)
        {
            if (node->left_vid == node->right_vid) { // both same partition
                unsigned partition = node->right_vid;
                FILE *file = pregraphFile[partition];
                fprintf(file, "NODE\t%"PRIu64"\t%hu\t%hu\n", ++(nodeCount[partition]),
                        (node->sequence.get_length()-g__FULLKMER_LENGTH+1), 0);
                node->sequence.PrintToFile(file);
                fputc('\n', file);
            } else { // different partitions
                {
                    unsigned partition = node->left_vid;
                    FILE *file = pregraphFile[partition];
                    fprintf(file, "NODE\t%"PRIu64"\t%hu\t%hu\n", ++(nodeCount[partition]),
                            (node->sequence.get_length()-g__FULLKMER_LENGTH+1), 1 /*, 1*/); // FIXME: boundary ID FIXME: arbitrarily chose left to emit
                    node->sequence.PrintToFile(file);
                    fputc('\n', file);
                }
                {
                    unsigned partition = node->right_vid;
                    FILE *file = pregraphFile[partition];
                    fprintf(file, "NODE\t%"PRIu64"\t%hu\t%hu\n", ++(nodeCount[partition]),
                            (node->sequence.get_length()-g__FULLKMER_LENGTH+1), 1 /*, 0*/); // FIXME: boundary ID FIXME: arbitrarily chose left to emit
                    node->sequence.PrintToFile(file);
                    fputc('\n', file);
                }
            }
        }
    };
    sg_for_each(sgraph, lambda2::emit_pregraph_sequence);


    //
    // read filtering
    //
    if (g__PGDIST_FILTER) {
        printf("Filter read set for %u partitions...\n", g__PGDIST_PARTITIONS); fflush(stdout);
        //  initialize sequence files
        for (unsigned i=0 ; i < g__PGDIST_PARTITIONS; ++i) {
            char newreads_filename[PATH_MAX+1];
            sprintf(newreads_filename, "%s/%u-Sequences", g__WORK_BASE_DIRECTORY, i);
            int newreads_filedes = open(newreads_filename, O_CREAT | O_TRUNC | O_WRONLY, 0666);
            readsFile[i] = fdopen(newreads_filedes, "w");
        }

        //
        // initialize the hash table
        //
        uint64_t buckets = (1 << g__PGDIST_FILTER);
        hashtable = (uint64_t*) calloc(buckets, sizeof(uint64_t));
        assert( hashtable != NULL );
        hashmask = buckets - 1;

        // TODO: alternatively, dump graph to disk and iterate over to fill read hashtable

        // fill hash table with partition mapping
        struct lambda3 {
            static void node_build_filter(SeqNode *node)
            {
                // fill hashtable
                lambda_fillhash fill;

                current_partition = node->left_vid;
                sequence_process_kmers(&node->sequence, g__FULLKMER_LENGTH, fill);
            
                // TODO: only process the kmers of a sequence once, by precomputing the partition bitmask   
                if (node->right_vid != node->left_vid) { 
                    current_partition = node->right_vid;
                    sequence_process_kmers(&node->sequence, g__FULLKMER_LENGTH, fill);
                }
            }
        };
        sg_for_each(sgraph, lambda3::node_build_filter);

        // filter reads
        char allreads_filename[PATH_MAX+1];
        sprintf(allreads_filename, "%s/Sequences", g__WORK_BASE_DIRECTORY);
        FILE *allreads = fopen(allreads_filename, "r");
        distribute_reads(g__PGDIST_PARTITIONS, g__FULLKMER_LENGTH, allreads);
        fclose(allreads);

        //  close sequence files
        for (unsigned i=0 ; i < g__PGDIST_PARTITIONS; ++i) {
            fclose(readsFile[i]);
        }
    }

    //
    // update pregraph node and read counts
    //
    for (unsigned i=0 ; i < g__PGDIST_PARTITIONS; ++i) {
        rewind(pregraphFile[i]);
        if (readCount[i] == 0) {
            fprintf(pregraphFile[i], "0x%08"PRIx64"", nodeCount[i]);
        } else {
            fprintf(pregraphFile[i], "0x%08"PRIx64"\t0x%08"PRIx64"", nodeCount[i], readCount[i]);
        }
        fclose(pregraphFile[i]);
    }
}

