#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@文件        :bam2rawgem.py
@说明        : trans bam file to raw gem file
@时间        :2022/03/18 22:49:50
@作者        :wjwei
@版本        :0.01
@邮箱        :wjwei9908@gmail.com
'''


import argparse
import logging
import os

import pysam
from tqdm import tqdm

# config the logging
logging.basicConfig(format='%(asctime)s - %(pathname)s[line:%(lineno)d] - %(levelname)s: %(message)s',
                    level=logging.DEBUG)

def getArgs():
    """receive the arguments

    Returns:
        args: args
    """
    # receive the arguments
    parser = argparse.ArgumentParser(description=__doc__, prog='gem_format.py')
    parser.add_argument('--gtf', action='store', dest='in_gtf', type=str, required=True,
                        help='your gtf file')
    parser.add_argument('--bam', action='store', dest='in_bam', type=str, required=True,
                        help='your gtf file')
    parser.add_argument('--indelsize', action='store', dest='indelsize', type=int, default=5,
                        help='allowd indel size to fix the hole')
    parser.add_argument('--gem', action='store', dest='out_gem', type=str, required=True,
                        help='output gem file name')
    return parser.parse_args()


def add_pos_flank(ivl, pos_str, size):
    """add flank to a interval list

    Args:
        ivl (list): interval list
        pos_str (str): postion
        size (int): flank size
    """
    pos = int(pos_str)
    for i in range(pos-size, pos+size+1):
        ivl.append(i)


def get_exon_ivl(gtf_file):
    """get exon interval list

    Args:
        gtf_file (str): gtf file path

    Returns:
        list: exon interval list
    """
    exon_ivl = []
    with open(args.in_gtf, 'r') as gtf_file_handler:
        for line in tqdm(gtf_file_handler, desc='Process gtf file'):
            if not line.startswith('#'):
                fields = line.rstrip().split('\t')
                chrom, feature_class, feature_type, start_str, end_str, junk, strand, junk, tags = fields
                if 'exon' in feature_type:
                    add_pos_flank(exon_ivl, start_str, 5)
                    add_pos_flank(exon_ivl, end_str, 5)
        # line = gtf_file_handler.readline()
        # while not line.startswith('#'):
        #     fields = line.rstrip().split('\t')
        #     chrom, feature_class, feature_type, start_str, end_str, junk, strand, junk, tags = fields
        #     if 'exon' in feature_type:
        #         add_pos_flank(exon_ivl, start_str, 5)
        #         add_pos_flank(exon_ivl, end_str, 5)
        #     line = gtf_file_handler.readline()

    return exon_ivl


def parse_cigar_tuple(cigartuples, pos, reverse, INDEL_SIZE):
    """parse cigar tuple of a read

    Args:
        cigartuples (tuples): cigar infomation
        pos (int): start position
        reverse (bool): if reverse
        INDEL_SIZE (int): allowed indel size for fix hole

    Returns:
        list: segments list
    """
    segments = []
    hole_to_remove = set()
    ref_skip = False
    clip5 = clip3 = 0
    p = pos
    cigartuples.reverse() if reverse else None
    for i, (operation_id, length) in enumerate(cigartuples):
        if operation_id == 0:  # vcy.CIGAR[operation_id] == "BAM_CMATCH"
            segments.append((p, p + length - 1))
            p += length
        # A splice || vcy.CIGAR[operation_id] == 'BAM_CREF_SKIP'
        elif operation_id == 3:
            ref_skip = True
            p += length
        # A deletion || cy.CIGAR[operation_id] == 'BAM_CDEL'
        elif operation_id == 2:
            if length <= INDEL_SIZE:
                try:
                    if cigartuples[i + 1][0] == 0 and cigartuples[i - 1][0] == 0:
                        hole_to_remove.add(len(segments) - 1)
                except IndexError:
                    pass
            p += length
        # bases at 5' or 3' are NOT part of the alignment || vcy.CIGAR[operation_id] == 'BAM_CSOFT_CLIP'
        elif operation_id == 4:
            if p == pos:
                clip5 = length  # At start of alignment
            else:
                # Must be at end of alignment vcy.CIGAR[operation_id] in ["BAM_CINS", "BAM_CHARD_CLIP"]
                clip3 = length
            p += length
        elif operation_id == 1:  # An insertion BAM_CINS
            if length <= INDEL_SIZE:
                try:
                    if cigartuples[i + 1][0] == 0 and cigartuples[i - 1][0] == 0:
                        hole_to_remove.add(len(segments) - 1)
                except IndexError:
                    pass
            # else do nothing
            # NOTE: maybe we should make so that the reads get discarded
        elif operation_id == 5:  # BAM_CHARD_CLIP
            print("Hard clip was encountered! All mapping are assumed soft clipped")

    # Merge segments separated by small insertions and deletions
    # NOTE maybe sorted is not required realy
    for a, b in enumerate(sorted(hole_to_remove)):
        segments[b - a] = (segments.pop(b - a)[0], segments[b - a][1])

    return segments


class ImportingError(Exception):
    """while file not existed or not a file or not readable, raise this!"""


class FileValidator:
    """validate the file """

    def validate(self, file_path: str):
        """check three posiable contions"""
        self._exists(file_path)
        self._is_file(file_path)
        self._is_readable(file_path)

    @classmethod
    def _exists(cls, file_path: str):
        """check if file exist"""
        if not os.path.exists(file_path):
            raise ImportingError(f"{file_path} does not exist")

    @classmethod
    def _is_file(cls, file_path: str):
        """check if is a file not a dir"""
        if not os.path.isfile(file_path):
            raise ImportingError(f"{file_path} is not a file")

    @classmethod
    def _is_readable(cls, file_path: str):
        """check if file readable"""
        if not os.access(file_path, os.R_OK):
            raise ImportingError(f"{file_path} is not readable")


def process_bamfile(bamfile, exon_ivl, indel_size, outgemfile):
    """process bamfile and output gem file

    Args:
        bamfile (str): bamfile path
        exon_ivl (list): flank added exon interval
        indel_size (int): allowed indel size to fix hole
        outgemfile (str): gemfile path

    Returns:
        4 statsis
    """
    multi_blocks = 0
    splice_beyond = 0
    intergeneic = 0
    unmapgene = 0
    bamhandler = pysam.AlignmentFile(bamfile, 'rb')
    with open(outgemfile, 'w') as out_gem_handler:
        for read in tqdm(bamhandler, desc='Process bamfile:'):
            status = 'UN'
            if read.is_unmapped or read.flag >= 256:
                continue
            else:
                refName = read.reference_name
                chrType = 'NORMAL' if refName.startswith('chr') else 'Contig'
                try:
                    gene = read.get_tag("GE:Z")
                    hit = read.get_tag("XF:Z")
                    spatial_pos = read.get_tag("CB:Z")
                    pos = read.reference_start + 1
                    segments = parse_cigar_tuple(
                        read.cigartuples, pos, reverse=read.is_reverse, INDEL_SIZE=indel_size)
                    if hit == 'INTRONIC':
                        status = 'Unspliced'
                    elif hit == 'EXONIC':
                        if 'N' in read.cigarstring:
                            starts, ends = zip(*segments)
                            if len(segments) >= 4:
                                multi_blocks += 1
                                continue
                            elif len(segments) == 3:
                                if starts[1] in exon_ivl and starts[2] in exon_ivl and ends[0] in exon_ivl and ends[1] in exon_ivl:
                                    status = 'Spliced'
                                else:
                                    splice_beyond += 1
                                    continue
                            else:
                                if starts[1] in exon_ivl and ends[0] in exon_ivl:
                                    status = 'Spliced'
                                else:
                                    splice_beyond += 1
                                    continue
                        else:
                            status = 'Ambiguous'
                    else:
                        intergeneic += 1
                        continue

                    writeline = '\t'.join(
                        [gene, spatial_pos, '1', status, chrType, str(segments)])
                    out_gem_handler.writelines(writeline + '\n')
                except:
                    unmapgene += 1
                    continue
    return multi_blocks, splice_beyond, intergeneic, unmapgene


if __name__ == '__main__':
    args = getArgs()
    # logging.info(
    #     "\033[32m Start to translate bem file {} to gem file {} \033[0m".format(
    #         args.in_bam,
    #         args.out_gem)
    #     )
    try:
        FileValidator().validate(args.in_gtf)
        logging.info(
            "\033[32m Start to read gtf file {} \033[0m".format(
                args.in_gtf,)
        )
        exon_ivl = get_exon_ivl(args.in_gtf)
    except ImportingError as error:
        raise ImportingError from error
    try:
        FileValidator().validate(args.in_bam)
        # FileValidator().validate(args.out_gem)
        logging.info(
            "\033[32m Start to parse bem file {} to gem file {} \033[0m".format(
                args.in_bam,
                args.out_gem)
        )
        multi_blocks, splice_beyond, intergeneic, unmapgene = process_bamfile(
            args.in_bam,
            exon_ivl,
            args.indelsize,
            args.out_gem)
        logging.info(
            "\033[32m Found: {} multi_blocks {} splice_beyond {} intergeneic  {} unmapgenes\033[0m".format(
                multi_blocks,
                splice_beyond,
                intergeneic,
                unmapgene)
        )
    except ImportingError as error:
        raise ImportingError from error
    logging.info(
        "\033[32m Al Done! Please check the gem file {} and format it \033[0m".format(
            args.out_gem)
    )

