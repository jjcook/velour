//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// mk_static.cpp
//
// partitioning method: load static mini-kmer assignment
//

#include "types.h"
#include "minikmer.h"

#include <math.h>
#include <set>

//
// version-dependent jazz
//

static int partfile_fd = 0;
static off_t partfile_length = 0;
static void * partfile_mmap = NULL;

Color computeKmerColor(Kmer fullKmer)
{
	return getMiniKmerColor( computeMiniKmer(fullKmer, g__FULLKMER_LENGTH, g__MINIKMER_LENGTH) );
}

void initializeMiniKmer(void)
{
	char exepath[PATH_MAX] = {0};
#ifdef __linux__
	readlink("/proc/self/exe", exepath, sizeof(exepath)); // exepath needs to be zero initialized
#else
    fprintf(stderr, "WARNING: Assumes static partition map in ./minikmer_ptables/\n");
	strcat(exepath, "./fake"); // FIXME: requires user to run exe from build directory
#endif
	char *basepath = dirname(exepath);

	char partfile[PATH_MAX] = {0};
	strncpy(partfile, basepath, PATH_MAX);
	sprintf(partfile, "%s/minikmer_ptables/%d.%d.mk", basepath, g__MINIKMER_LENGTH, g__PARTITION_COUNT);

   	partfile_fd = open(partfile, O_RDONLY);
	if (partfile_fd == -1) {
		fprintf(stderr,"failed to open file: %s\n", partfile); 
		exit(EXIT_FAILURE);
	}

	struct stat file_stat;
	if (fstat(partfile_fd, &file_stat) != 0) {
		fprintf(stderr,"failed to stat file: %s\n", partfile); 
		exit(EXIT_FAILURE);
	}

	partfile_length = file_stat.st_size;
	partfile_mmap = mmap(0, partfile_length, PROT_READ, MAP_PRIVATE, partfile_fd, 0);
	if (partfile_mmap == (void *)-1) {              
		fprintf(stderr, "failed to mmap file: %s\n", partfile); 
		exit(EXIT_FAILURE);
	}

	printf("Mapped minikmer mapping file: %s\n", partfile);

	miniKmerGraph = static_cast<struct minikmer_node_s *>(static_cast<void *>((static_cast<char*>(partfile_mmap) + sizeof(uint64_t) + sizeof(uint64_t))));

  	MINIKMER_GRAPH_SIZE = 1UL << (2 * g__MINIKMER_LENGTH);
	//assert( partfile_length - sizeof(uint64_t) - sizeof(uint64_t) == MINIKMER_GRAPH_SIZE * sizeof(struct minikmer_node_s) );

	COLORMIN = 0;
	COLORMAX = COLORMIN + g__PARTITION_COUNT - 1;

	printf("MKG Initialization:  mini-k = %d  full-k = %d  colors = %u  totalNodes = %" PRIuPTR"\n",
	 	g__MINIKMER_LENGTH, g__FULLKMER_LENGTH, g__PARTITION_COUNT, MINIKMER_GRAPH_SIZE );
}

void destroyMiniKmer(void)
{
	if (partfile_fd != 0) {
		munmap(partfile_mmap, partfile_length);
		close(partfile_fd);
		partfile_fd = 0;
	}
}

void printMiniKmerStats(void)
{
	printMiniKmerGraphAdjacencyStats();
}

