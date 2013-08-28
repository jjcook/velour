=============================================================================
VELOUR IMPLEMENTATION DETAILS AND DISCLAIMERS
=============================================================================

=== DISCLAIMERS ===

- I have not had the time to refactor many of the graph inspection & traversal
  mnemonics, nor remove deprecated code that is commented out nor ifdef'd out;
  overall, there is significant opportunity for code clean up.

- Multi-threaded execution is experimental.  Specifically, it has not been
  fully tested (nor validated with Thread Checker) to ensure data-race
  freedom.

=== IMPLEMENTATION DETAILS ===

<< todo >>

- A back-end to apply Velvet's correctGraph algorithm to partitions of the
  PreGraph is mostly (but not fully) implemented; this process is described in
  the appendix of my disseration, called localization.

=== FUTURE FEATURE IDEAS ===

- large k-mer support (greater than 31)
- large sequence node support, currently limited to 65535 bases

=== FUTURE OPTIMIZATION IDEAS ===

- for bubble popping, use a fast banded comparison algorithm


NOTE: below are some old todo notes, some may have already been implemented, others may no longer apply or make sense

- paralleizing redistribution over nfs: use netcat command 'nc' on named fifo pipes on each side of fifo for each node

!! only need to error correct / concat the component that we just loaded, not the entire graph!!!

- fflush stdout and stderr before printing error messages
- hash table size function of input/genome size
- global variables: FULLKMER_MASK and MINIKMER_MASK
- check results vs LOOM to check for bugs!
- signal handler that flushes stdout/stderr -- or use setvbuf? detect if terminal?
- madvise, fadvise, mmap(MAP_POPULATE) for smaller files, ... POSIX_FADV_SEQUENTIAL | POSIX_FADV_NOREUSE
- input-agnostic tuning script that generates a machine/mount-specific config file for velour
- get/setrlimit for resource limits
- malloc/bzero -> calloc
- read() vs mmap
- cleanup: convert divide/mod by 4 into logical shift/and?
- cleanup: test + error message after every allocation?
- use typeof operator when doing sizeof? err no, that's not how it works
- sizeof(variable) instead of sizeof(type)
- NORANDOM option, so no funny error messages about /dev/urandom
- fread/fwrite return values
- parallel seq graph insertion scheme: insert in normal hashtable as done for kmer graph, if actually inserted, then insert into rc hashtable
    (do I need to canonicalize the sequence?  in this scheme, YES)
- mmap: osx/bsd: MAP_NOCACHE
- flush after printing to terminal --> disable buffering

- kg/sg_node color fields: size to needed # of partitions, jam flags into leftover space if available
- distribution function: return true if need another pass over the read set

- dump irreducible subcomponents to final bucket?
- newer TBB has better generic atomics support now
- "verify" functions: check only happening when supposed to
- "explicit" class constructors
- check if sequence gets longer than 'unsigned' bases
- disable exceptions? and RTTI
- use TBB allocator: convert std::vector<T> to std::vector<T,tbb::scalable_allocator<T> >

- prefer appends to prepends somehow
- posix_fallocate() to pre-allocate file space -- if using filesytem with 'extents' which auto-magically happens to zero it too
  -- fallocate() is syscall, used if available by posix_fallocate()

- GNU:  __attribute__((always_inline)) 

- SIGABRT / SIGSEGV handler to invoke debugger (environment variable enable this) (will this work if stdout is redirected!?)
  - example:
        void sigabrt_handler(int signum) {
            exit(EXIT_FAILURE);
        }

        void cleanup(void) {
            /* delete temporary files, restore consistent state, etc */
        }

        int main(void) {
            atexit(cleanup);
            signal(SIGABRT, sigabrt_handler);
            /* ... */
            assert(/* something bad didn't happen */);
            /* ... */
        }

- makefile: check TBB environment variables so compilation doesn't vomit
- open(): O_SYMLINK follow symlinks
- change ATOMIZE to ATOMIC_LOAD / ATOMIC_STORE for clarity where appropriate
- memory management: do not garbage collect graph constantly if using swap file, as that will hurt?
- callorOrExit / mallocOrExit

- use large page for mini-kmer lookup table for hash table localization
- increase size of partition file buffer, as looming can be easily disk limited
- simple statistic of size (nodes) of largest component observed during flowing bucket
- catch case of non-existant inbox???
- peak live memory in allocator?
