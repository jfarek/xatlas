#include "Bam.hpp"
#include <cstring>
#include <cstdlib>

/**
 * Handle BAM or CRAM file
 */

Bam::Bam()
    : _idx(nullptr),
      _thread_pool(nullptr),
      _sf(nullptr),
      _iter(nullptr),
      _hdr(nullptr),
      _status(BAM_OKAY)
{
}

Bam::Bam(const char *sf_fn, const char *ref_fn, uint8_t nthreads)
    : _idx(nullptr),
      _thread_pool(nullptr),
      _sf(nullptr),
      _iter(nullptr),
      _hdr(nullptr),
      _status(BAM_OKAY)
{
    size_t fn_strlen = strlen(sf_fn);
    size_t ref_fn_strlen = strlen(ref_fn);
    int set_fai_status;
    _status = BAM_OKAY;

    if (strcmp(sf_fn + fn_strlen - 3, "bam") == 0) {
        // load BAM file
        _sf = sam_open(sf_fn, "rb");
    } else if (strcmp(sf_fn + fn_strlen - 4, "cram") == 0) {
        // load CRAM file
        _sf = sam_open(sf_fn, "rc");

        // load reference index for reading cram
        char *fai_fn = new char[ref_fn_strlen + 5];
        std::strncpy(fai_fn, ref_fn, ref_fn_strlen);
        fai_fn[ref_fn_strlen] = '\0';
        std::strncat(fai_fn, ".fai", 4);

        set_fai_status = hts_set_fai_filename(_sf, fai_fn);
        delete[] fai_fn;
        if (set_fai_status != 0) {
            _status = BAM_BAD_REF_INDEX;
            return;
        }
    } else {
        // load as BAM format
        _sf = sam_open(sf_fn, "rb");
    }

    if (_sf == nullptr) {
        _status = BAM_BAD_FILE;
    } else if ((_hdr = sam_hdr_read(_sf)) == nullptr) {
        _status = BAM_BAD_HEADER;
    } else if ((_idx = sam_index_load(_sf, sf_fn)) == nullptr) {
        _status = BAM_BAD_INDEX;
    } else if (nthreads > 1) {
        // activate htslib i/o multithreading
        _thread_pool = (htsThreadPool *)std::malloc(sizeof(htsThreadPool));
        if (_thread_pool != nullptr) {
            _thread_pool->qsize = 0;
            _thread_pool->pool = hts_tpool_init(nthreads);
            hts_set_opt(_sf, HTS_OPT_THREAD_POOL, _thread_pool);
        }
    }
}

Bam::Bam(const Bam &) = default;

Bam::~Bam()
{
    hts_itr_destroy(_iter);
    hts_idx_destroy(_idx);
    bam_hdr_destroy(_hdr);
    sam_close(_sf);
    if (_thread_pool != nullptr) {
        hts_tpool_destroy(_thread_pool->pool);
        std::free(_thread_pool);
    }
}

void Bam::set_iter(const char *region_str)
{
    hts_itr_destroy(_iter);
    _iter = sam_itr_querys(_idx, _hdr, region_str);
}
