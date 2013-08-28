//
// This file is part of Velour and is distributed under the University of
// Illinois Open Source License.  See LICENSE.txt for details.
//

//
// allocators.cpp
//

#include "types.h"

#ifdef VELOUR_TBB
#  ifndef USE_TBB_ALLOC
#    include <tbb/spin_mutex.h> // which TBB has already
#  endif
#  include <tbb/atomic.h> // Used for tracking wasted space in the allocator
#endif

//********************************************************************************
//*************************    Large Page Support    *****************************
//********************************************************************************

#ifndef MAP_ANONYMOUS
#define MAP_ANONYMOUS MAP_ANON
#endif

// TODO: page size should be determined not assumed?
// Use getpagesize() from unistd.h
#define PAGE_SIZE (4*1024)
#define HUGE_PAGE_SIZE (2*1024*1024)
// TODO: allocation size should _not_ be fixed!?
#define ALLOC_SIZE (32UL*1024*1024)

typedef struct allocator {
  int largepage_fd;  
  char *start;
  char *current;
  char *end;
  char *filename;
  struct allocator *next;
} allocator_t;

static allocator_t g_allocator;

#ifdef VELOUR_TBB

#   define TLS_ALLOC_SIZE (1UL*1024*1024) // Size of per-thread pools

    typedef struct thread_allocator {
        char *start;
        char *current;
        char *end;
    } thread_allocator_t;

    static __thread thread_allocator_t tls_allocator; // One of these exists per-thread

#endif // VELOUR_TBB
    
#ifdef VELOUR_TBB    
static tbb::atomic<size_t> WASTED_SPACE; // Keep track of the amount of space the TLS allocator wastes
#else
static size_t WASTED_SPACE;
#endif

void
velour_alloc_init() {
#ifdef USE_TBB_ALLOC
  // scalable_malloc requires no initialization
  return;
#else

  WASTED_SPACE = 0;
  
  g_allocator.start = 0;
  g_allocator.next = NULL;
  g_allocator.filename = NULL;

  char *filename = (char *)malloc((PATH_MAX+1) * sizeof(char));
  strncpy(filename, "/mnt/hugetlbfs/Velvet.recyclebin.", PATH_MAX);
  pid_t pid = getpid();
  sprintf(&filename[strlen(filename)], "%d", pid);
  //printf("%s\n", filename);
  g_allocator.largepage_fd = open(filename, O_CREAT | O_RDWR, 0755);

  if (g_allocator.largepage_fd != -1) { 
    g_allocator.start = (char *)mmap(0x0, ALLOC_SIZE, PROT_READ | PROT_WRITE, 
                                     MAP_SHARED, g_allocator.largepage_fd, 0);
    if (g_allocator.start == (char *)-1) {
      printf("falling back on standard allocator\n");
      g_allocator.start = 0; // fall back on the standard allocator
    } else {
      g_allocator.filename = filename;
      printf("using large pages....\n");
    }
  }
  if (g_allocator.start == 0) { // fall back on the default allocator
    g_allocator.start = (char *)mmap(0, ALLOC_SIZE, PROT_READ | PROT_WRITE, 
                                     MAP_SHARED|MAP_ANONYMOUS, 0, 0);
    // g_allocator.start = (char *)malloc(ALLOC_SIZE);
    if ((g_allocator.start == NULL) || g_allocator.start == (char *)MAP_FAILED) {
      printf("failed memory allocation for hash table and pre-nodes\n");
      exit(1);
    }
  }
  g_allocator.current = g_allocator.start;
  g_allocator.end =     g_allocator.start + ALLOC_SIZE;
  if (g_allocator.filename == NULL) {
	 free(filename);
  }
#endif
}

void
velour_alloc_done() {
#ifdef USE_TBB_ALLOC
  // TODO: Find a way to do clean free-ing of memory from scalable_alloc
  return;
#else

  for (allocator_t *trav = &g_allocator ; trav != NULL ; /* nothing */) {
	 munmap((void *)trav->start, (int)(trav->end - trav->start));
	 if (trav->largepage_fd != -1) {
		close(trav->largepage_fd);
	 }
	 if (trav->filename != NULL) {
		unlink(trav->filename);
		// free(trav->filename);
	 }
	 allocator_t *temp = trav;
	 trav = trav->next;
	 if (temp != &g_allocator) {
		free(temp);
	 }
  }
  g_allocator.next = NULL;
  g_allocator.largepage_fd = -1;

  printf("Allocator: %"PRIuPTR" bytes wasted\n", reinterpret_cast<const uintptr_t&>(WASTED_SPACE)); // or "+0" to get correct value?
#endif
}

static void
velour_alloc_reinit() {
  allocator_t *temp = (allocator_t *)malloc(sizeof(allocator_t));
  memcpy(temp, &g_allocator, sizeof(allocator_t));
  g_allocator.next = temp;
  g_allocator.largepage_fd = -1;
  g_allocator.filename = NULL;
  char *addr = (char *)mmap(temp->end, ALLOC_SIZE/2, PROT_READ | PROT_WRITE, 
                            MAP_SHARED|MAP_ANONYMOUS, 0, 0);
  if (addr == (char *)MAP_FAILED) {
	 printf("allocator ran out of memory\n");
	 exit(1);
  }
  g_allocator.start =   addr;
  g_allocator.current = addr;
  g_allocator.end =     addr + (ALLOC_SIZE/2);
}

#ifdef VELOUR_TBB
#  ifndef USE_TBB_ALLOC
// Big spinlock for the bulk, centralized allocator
static tbb::spin_mutex the_alloc_lock;
#  endif
#endif

void *
velour_alloc(size_t sz) {
#ifdef USE_LIBC_ALLOC
  return malloc(sz);
#else

#ifdef USE_TBB_ALLOC
  // Skip all the bookkeeping, just use TBB's alloc
  return scalable_malloc(sz);
#else

#ifdef VELOUR_TBB // note:  didn't get here unless we're using velour_alloc
  tbb::spin_mutex::scoped_lock lock;
  lock.acquire(the_alloc_lock);
#endif
  if ((g_allocator.current + sz) > g_allocator.end) {
	 velour_alloc_reinit();
  }
  void *ret_val = (void *)g_allocator.current;
  g_allocator.current += sz;
#ifdef VELOUR_TBB
  lock.release();
#endif
  return ret_val;
#endif // not TBB_ALLOC
#endif // not LIBC_ALLOC
}

void *
velour_realloc(void * old_ptr, size_t old_sz, size_t new_sz) {
#ifdef USE_LIBC_ALLOC
  return realloc(old_ptr, new_sz);
#else

#ifdef USE_TBB_ALLOC
  // Skip all the bookkeeping, just use TBB's alloc
  return scalable_realloc(old_ptr, new_sz);
#else

#ifdef VELOUR_TBB // note:  didn't get here unless we're using velour_alloc
  tbb::spin_mutex::scoped_lock lock;
  lock.acquire(the_alloc_lock);
#endif
    void *ret_val = NULL;
    if ((g_allocator.current - old_sz) == old_ptr) {
        if (old_sz > new_sz) {
            g_allocator.current -= (old_sz - new_sz);
            ret_val = old_ptr;
        } else {
            if ((g_allocator.current + (new_sz - old_sz)) > g_allocator.end) {
                velour_alloc_reinit();
                ret_val = (void *)g_allocator.current;
                memcpy(ret_val, old_ptr, old_sz);
                g_allocator.current += new_sz;
            } else {
                g_allocator.current += (new_sz - old_sz);
                ret_val = old_ptr;
            }
        }
    } else {
        if ((g_allocator.current + new_sz) > g_allocator.end) {
            velour_alloc_reinit();
        }
        ret_val = (void *)g_allocator.current;
        memcpy(ret_val, old_ptr, old_sz);
        g_allocator.current += new_sz;
    }
#ifdef VELOUR_TBB
  lock.release();
#endif
  return ret_val;
#endif // not TBB_ALLOC
#endif // not LIBC_ALLOC
}

// this could be improved
void * velour_calloc(size_t sz)
{
#ifdef USE_LIBC_ALLOC
  return calloc(1, sz);
#else

#ifdef USE_TBB_ALLOC
  // Skip all the bookkeeping, just use TBB's alloc
  return scalable_calloc(1, sz);
#else

    void *mem = velour_alloc(sz);

    uintptr_t uintptr_limit = sz / sizeof(uintptr_t);
    for (size_t i=0; i < uintptr_limit; ++i) {
        static_cast<uintptr_t *>(mem)[i] = 0;
    }

    size_t mod = sz % sizeof(uintptr_t);
    if (mod > 0) {
        char *offset = static_cast<char*>(mem) + (uintptr_limit * sizeof(uintptr_t));
        char *limit = offset + mod;
        for (; offset < limit; ++offset) {
            *offset = 0;
        }
    }

    return mem;
#endif // not TBB_ALLOC
#endif // not LIBC_ALLOC
}

#ifdef VELOUR_TBB
// Thread-local allocation; tls_allocator is unique per-thread
// This allocator is built on top of the velour_alloc
void *
tls_velour_alloc(size_t sz) {
#ifdef USE_LIBC_ALLOC
  return malloc(sz);
#else

#ifdef USE_TBB_ALLOC
  // Skip all the bookkeeping, just use TBB's alloc
  return scalable_malloc(sz);
#else
  // check if we can service this request at all
  if(sz > TLS_ALLOC_SIZE){ 
    //nope, kick it up the chain
    return velour_alloc(sz);
  }
  // Figure how much space remains - this works for initialization as well
  size_t remaining_bytes = tls_allocator.end - tls_allocator.current;
  if(remaining_bytes < sz) {
    // refill the local buffer
    WASTED_SPACE += remaining_bytes;
    tls_allocator.current = tls_allocator.start = (char *)velour_alloc(TLS_ALLOC_SIZE);
    tls_allocator.end = tls_allocator.start + TLS_ALLOC_SIZE;
  }
  void *ret_val = (void *)tls_allocator.current;
  tls_allocator.current += sz;
  return ret_val;
#endif // not TBB_ALLOC
#endif // not LIBC_ALLOC
}
#endif // VELOUR_TBB

//********************************************************************************
//**************************    "Pool" Allocation    *****************************
//********************************************************************************

static unsigned g_num_allocator_pools = 0; // TODO unused?
static unsigned g_allocator_pool_shift = 0; // TODO unused?
static allocator_t *g_allocator_pools = NULL;

void
velour_allocate_pools() {
  size_t pages_avail = (g_allocator.end - g_allocator.current)/HUGE_PAGE_SIZE;
  // find largest power of two number of pools, each 1 page large
  unsigned ntemp = 0;
  while (pages_avail >= (unsigned)(2 << ntemp)) {
	 ntemp ++;
  }
  g_num_allocator_pools = 1 << ntemp;
  g_allocator_pool_shift = HASH_BITS - ntemp;
  g_allocator_pools = (allocator_t *)malloc(g_num_allocator_pools * sizeof(allocator_t));
  for (unsigned i = 0 ; i < g_num_allocator_pools ; i ++) {
	 char *temp = g_allocator.current + (i*HUGE_PAGE_SIZE);
	 g_allocator_pools[i].start = temp;
	 g_allocator_pools[i].current = temp;
	 g_allocator_pools[i].end = temp + HUGE_PAGE_SIZE;
  }
  g_allocator.current += g_num_allocator_pools * HUGE_PAGE_SIZE;
}

static void 
velour_refill_pool(unsigned index) {
  assert(index < g_num_allocator_pools);
  off_t length =  g_allocator_pools[index].end - g_allocator_pools[index].start;
  length = (length > PAGE_SIZE) ? length/2 : PAGE_SIZE;
  g_allocator_pools[index].start   = (char *) velour_alloc(length);
  g_allocator_pools[index].current = g_allocator_pools[index].start;
  g_allocator_pools[index].end     = g_allocator_pools[index].start + length;
}

void *
velour_pool_alloc(size_t sz, uintptr_t hash) {
  unsigned index = hash >> g_allocator_pool_shift;
  if ((g_allocator_pools[index].current + sz) > g_allocator_pools[index].end) {
	 velour_refill_pool(index);
  }
  void *ret_val = (void *)g_allocator_pools[index].current;
  g_allocator_pools[index].current += sz;
  return ret_val;
}

