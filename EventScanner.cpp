#include <algorithm>
#include <deque>
#include "EventScanner.hpp"

EventScanner::EventScanner()
{
}

EventScanner::EventScanner(int32_t near_end, double snp_max_sub, double snp_max_gap)
    : _near_end(near_end),
      _snp_max_sub(snp_max_sub),
      _snp_max_gap(snp_max_gap),
      _last_call_indel(-1),
      _last_call_snp(-1)
{
}

EventScanner::EventScanner(const EventScanner &) = default;

EventScanner::~EventScanner() = default;

void EventScanner::reset()
{
    _indels.clear();
    _snps.clear();
    _last_call_indel = -1;
    _last_call_snp = -1;
}

int32_t EventScanner::left_align(const ReferenceSequence &refseq, std::string& seq, int32_t indel_pos, int32_t read_start_pos)
{
    int32_t align_pos = indel_pos;
    size_t len = seq.length();
    size_t len_1 = len - 1;
    uint32_t ref_idx = (uint32_t)indel_pos - 1;
    uint32_t seq_idx = 0;

    while (align_pos > read_start_pos && refseq._seq[ref_idx] == seq[len_1 - seq_idx]) {
        --align_pos;
        --ref_idx;
        seq_idx = (seq_idx + 1) % len;
    }

    return align_pos;
}

void EventScanner::collect_indels(
        bam1_t *read,
        const ReferenceSequence &refseq,
        CoverageCounter &coverages,
        bed_coord_t &segment)
{
    bool strand, first_op;
    uint8_t *bseq;
    uint32_t var_count, c_op, c_type, *cigar, *end_cigar;
    int32_t c_len, pos1, pos2, end_pos, indel_pos, left_align_pos;
    int32_t q_pos, r_pos, leading_sclip, tailing_sclip;
    std::string seq;
    AReadsIndelList this_reads_indels;

    first_op = true;
    var_count = 0;
    cigar = bam_get_cigar(read);
    end_cigar = cigar + read->core.n_cigar - 1;
    indel_pos = read->core.pos;
    q_pos = 0;
    r_pos = 0;
    leading_sclip = (bam_cigar_op(*cigar) == BAM_CSOFT_CLIP)
                        ? bam_cigar_oplen(*cigar)
                        : 0;
    tailing_sclip = (bam_cigar_op(*end_cigar) == BAM_CSOFT_CLIP)
                        ? bam_cigar_oplen(*end_cigar)
                        : 0;
    bseq = bam_get_seq(read);
    strand = (bam_is_rev(read) == 0);

    // for adjust reference coverage, deque of (pos, size)
    std::deque<std::pair<int32_t, int32_t>> coverage_rr_vec, coverage_del_vec;

    for (; cigar <= end_cigar; ++cigar) {
        c_op = bam_cigar_op(*cigar);
        c_len = (int32_t)bam_cigar_oplen(*cigar);

        // coverage and events
        switch (c_op) {
        case BAM_CMATCH:
            //coverages.add_coverage_range(indel_pos, COVSIDX_INDEL_RR, c_len);
            coverage_rr_vec.push_back(std::make_pair(indel_pos, c_len));
            coverages.add_coverage_range(indel_pos, COVSIDX_INDEL_DP, c_len);
            pos1 = q_pos + leading_sclip;
            pos2 = read->core.pos + r_pos;
            end_pos = pos1 + c_len;
            while (pos1 < end_pos) {
                if (seq_nt16_str[bam_seqi(bseq, pos1)] != refseq._seq[pos2]) {
                    ++var_count;
                }
                ++pos1;
                ++pos2;
            }
            break;
        case BAM_CREF_SKIP:
        case BAM_CEQUAL:
            //coverages.add_coverage_range(indel_pos, COVSIDX_INDEL_RR, c_len);
            coverage_rr_vec.push_back(std::make_pair(indel_pos, c_len));
            coverages.add_coverage_range(indel_pos, COVSIDX_INDEL_DP, c_len);
            break;
        case BAM_CDIFF:
            //coverages.add_coverage_range(indel_pos, COVSIDX_INDEL_RR, c_len);
            coverage_rr_vec.push_back(std::make_pair(indel_pos, c_len));
            coverages.add_coverage_range(indel_pos, COVSIDX_INDEL_DP, c_len);
            var_count += c_len;
            break;
        case BAM_CINS:
            pos1 = q_pos + leading_sclip;
            end_pos = pos1 + c_len;
            seq.clear();
            while (pos1 < end_pos) {
                seq.push_back(seq_nt16_str[bam_seqi(bseq, pos1)]);
                ++pos1;
            }
            left_align_pos = left_align(refseq, seq, indel_pos, read->core.pos);

            // rotate seq if different left align position
            if (left_align_pos != indel_pos) {
                std::rotate(seq.begin(), seq.begin() + c_len - (indel_pos - left_align_pos) % c_len, seq.end());
            }

            if (left_align_pos >= segment.first &&
                left_align_pos <= segment.second)
            {
                this_reads_indels.push_back(
                    IndelEvent(
                        left_align_pos,
                        c_len,
                        q_pos,
                        seq,
                        strand,
                        read->core.qual));
            }
            if (!first_op && indel_pos >= 1) {
                coverages.add_coverage(left_align_pos - 1, COVSIDX_INDEL_RR_INS);
            }
            ++var_count;
            break;
        case BAM_CDEL:
            pos1 = r_pos + leading_sclip;
            end_pos = pos1 + c_len;
            seq.clear();
            while (pos1 < end_pos) {
                seq.push_back(refseq._seq[pos1]);
                ++pos1;
            }
            left_align_pos = left_align(refseq, seq, indel_pos, read->core.pos);
            seq.clear();
            if (left_align_pos >= segment.first &&
                left_align_pos <= segment.second &&
                left_align_pos + c_len - 1 <= segment.second)
            {
                this_reads_indels.push_back(
                    IndelEvent(
                        left_align_pos,
                        c_len,
                        q_pos,
                        seq,
                        strand,
                        read->core.qual));
            }
            coverages.add_coverage_range(indel_pos, COVSIDX_INDEL_DP, c_len);
            coverage_del_vec.push_back(std::make_pair(left_align_pos, c_len));
            ++var_count;
            break;
        default:
            break;
        }
        first_op = false;

        c_type = bam_cigar_type(*cigar);
        if ((c_type & 0x01) != 0 && c_op != BAM_CSOFT_CLIP) {
            q_pos += c_len;
        }
        if ((c_type & 0x02) != 0) {
            r_pos += c_len;
            indel_pos += c_len;
        }
    }

    // add RR coverage, adjusting for deletions
    for (; !coverage_rr_vec.empty(); coverage_rr_vec.pop_front())
    {
        auto &rr_itr = coverage_rr_vec.front();
        int32_t next_del_pos, next_del_size;
        int32_t pos = rr_itr.first;
        int32_t end = pos + rr_itr.second;

        for (; !coverage_del_vec.empty(); coverage_del_vec.pop_front()) {
            auto &del_itr = coverage_del_vec.front();
            next_del_pos = del_itr.first;
            if (next_del_pos >= end) {
                break;
            }
            coverages.add_coverage_range(pos, COVSIDX_INDEL_RR, next_del_pos - pos);
            next_del_size = del_itr.second;
            pos += next_del_size;
            end += next_del_size;
        }

        coverages.add_coverage_range(pos, COVSIDX_INDEL_RR, end - pos);
    }

    if (!this_reads_indels.empty()) {
        uint8_t *bqual = bam_get_qual(read) + leading_sclip;
        int32_t rlen = read->core.l_qseq - leading_sclip - tailing_sclip;
        double var_rate_gap_and_mismatch = (double)var_count / rlen;
        qual_t qual_sum;

        for (auto &iv_it : this_reads_indels) {
            end_pos = iv_it._q_pos;
            if (!iv_it._seq.empty()) {
                end_pos += iv_it._var_len;
            }

            pos1 = iv_it._q_pos - _near_end - 2;
            if (pos1 < 0) {
                pos1 = 0;
            }

            pos2 = end_pos + _near_end;
            if (pos2 >= rlen) {
                pos2 = rlen - 1;
            }

            qual_sum = 0;
            for (q_pos = pos1; q_pos <= pos2; ++q_pos) {
                qual_sum += (qual_t)bqual[q_pos];
            }

            iv_it._near_read_end_count = (iv_it._q_pos + 1 <= _near_end || rlen - end_pos < _near_end + 2 ? 1 : 0);
            iv_it._avg_nbq = qual_sum / (pos2 - pos1 + 1);
            iv_it._var_rate_gap_and_mismatch = var_rate_gap_and_mismatch;

            const auto &iv_search = _indels.find(iv_it._var_start);
            if (iv_search == _indels.end()) {
                _indels[iv_it._var_start].insert(std::make_pair(iv_it._id, iv_it));
            } else {
                const auto &ivg_it = iv_search->second.find(iv_it._id);
                if (ivg_it == iv_search->second.end()) {
                    iv_search->second.insert(std::make_pair(iv_it._id, iv_it));
                } else {
                    ivg_it->second.add_indel_event(iv_it);
                }
            }
        }
    }
}

/**
 * Sum of base quality score within range [p-5, p+4]
 */
inline double EventScanner::snp_nqs(int32_t p, int32_t qlen, const uint8_t *bqual)
{
    int32_t start, end;
    qual_t qual_sum;

    start = (p < 5) ? 0 : p - 5;
    end = (p > qlen - 4) ? qlen : p + 4;
    qual_sum = 0;

    for (int32_t qpos = start; qpos <= end; ++qpos) {
        qual_sum += (qual_t)bqual[qpos];
    }

    return (double)qual_sum / 10.0; // NOTE adjust for near read ends
}

void EventScanner::collect_snps(bam1_t *read, const ReferenceSequence &refseq, CoverageCounter &coverages, bed_coord_t &segment)
{
    bool strand;
    char refbase, allele;
    uint8_t snp_nt_gap, nm_minus_gap, snp_nt_sub, snp_nt_snp, *bseq, *bqual;
    uint32_t *cigar, *cigar0, *end_cigar;
    int32_t c_len, len, pos1, pos2, r_pos, q_pos;
    AReadsSnpList this_reads_snps;

    snp_nt_gap = 0;
    len = read->core.l_qseq;
    cigar0 = bam_get_cigar(read);
    cigar = cigar0;
    end_cigar = cigar + read->core.n_cigar - 1;

    // snp sub and gap
    while (cigar <= end_cigar) {
        switch (bam_cigar_op(*cigar)) {
        case BAM_CINS:
        case BAM_CDEL:
            snp_nt_gap += bam_cigar_oplen(*cigar);
            break;
        default:
            break;
        }

        ++cigar;
    }

    nm_minus_gap = bam_aux2i(bam_aux_get(read, "NM")) - snp_nt_gap;
    snp_nt_sub = nm_minus_gap;

    if ((double)snp_nt_sub / len > _snp_max_sub ||
        (double)snp_nt_gap / len > _snp_max_gap)
    {
        return;
    }

    // coverage and events
    strand = (bam_is_rev(read) == 0);
    r_pos = read->core.pos;
    q_pos = 0;
    snp_nt_snp = 0;
    bseq = bam_get_seq(read);
    bqual = bam_get_qual(read);

    // collect snps from cigar M mismatches
    for (cigar = cigar0; cigar <= end_cigar; ++cigar) {
        switch (bam_cigar_op(*cigar)) {
        case BAM_CMATCH:
            c_len = bam_cigar_oplen(*cigar);
            if (nm_minus_gap == 0) {
                // no mismatches remaining in read, only need to record coverage
                coverages.add_coverage_range(r_pos, COVSIDX_SNP_RR, c_len);
                coverages.add_coverage_range(r_pos, COVSIDX_SNP_DP, c_len);
                r_pos += c_len;
                q_pos += c_len;
            } else {
                // find snps
                pos1 = r_pos;
                len = r_pos + c_len;

                while (r_pos < len) {
                    refbase = refseq._seq[r_pos];
                    allele = seq_nt16_str[bam_seqi(bseq, q_pos)];

                    // snp found
                    if (refbase != allele) {
                        ++snp_nt_snp;
                        --nm_minus_gap;

                        // coverage for bases before snp
                        if (r_pos > pos1) {
                            pos2 = r_pos - pos1;
                            coverages.add_coverage_range(pos1, COVSIDX_SNP_RR, pos2);
                            coverages.add_coverage_range(pos1, COVSIDX_SNP_DP, pos2);
                        }
                        pos1 = r_pos + 1;

                        // add snp
                        if (!(refbase == 'N' || allele == 'N') &&
                            r_pos >= segment.first &&
                            r_pos <= segment.second)
                        {
                            double relpos = (double)q_pos / (double)read->core.l_qseq;
                            this_reads_snps.push_back(
                                std::make_pair(
                                    r_pos,
                                    SnpEvent(
                                        allele,
                                        strand,
                                        (qual_t)bqual[q_pos],
                                        relpos,
                                        this->snp_nqs(q_pos, read->core.l_qseq, bqual))));
                        }
                    }

                    ++r_pos;
                    ++q_pos;
                }

                // coverage for remaining bases
                if (r_pos > pos1) {
                    pos2 = r_pos - pos1;
                    coverages.add_coverage_range(pos1, COVSIDX_SNP_RR, pos2);
                    coverages.add_coverage_range(pos1, COVSIDX_SNP_DP, pos2);
                }
            }
            break;
        case BAM_CINS:
        case BAM_CSOFT_CLIP:
            q_pos += bam_cigar_oplen(*cigar);
            break;
        case BAM_CDEL:
            c_len = bam_cigar_oplen(*cigar);
            coverages.add_coverage_range(r_pos, COVSIDX_SNP_DP, c_len);
            r_pos += c_len;
            break;
        /* probably should process these too
        case BAM_CDIFF:
        case BAM_CEQUAL:
        case BAM_CPAD:
        case BAM_CREF_SKIP:
        case BAM_CHARD_CLIP:
            break;
        */
        default:
            break;
        }
    }

    // add each snp
    if (snp_nt_sub == snp_nt_snp) {
        for (const auto &sv : this_reads_snps) {
            const auto &sv_it = _snps.find(sv.first);
            if (sv_it == _snps.end()) {
                _snps[sv.first].insert(std::make_pair(sv.second._allele_base, sv.second));
            } else {
                const auto &svg_it = sv_it->second.find(sv.second._allele_base);
                if (svg_it == sv_it->second.end()) {
                    sv_it->second.insert(std::make_pair(sv.second._allele_base, sv.second));
                } else {
                    svg_it->second.add_snp_event(sv.second);
                }
            }
        }
    }
}
