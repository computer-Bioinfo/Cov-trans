
import pysam
import pandas as pd
import sys
import argparse
from collections import Counter
import gffutils



class SpliceReads():
    def __init__(self, reference: str,  contig=None):
        self.CIGAR = ["M", "I", "D", "N", "S", "H", "P", "=", "X", "B"]
        self.sgRNAs = {}
        self.ORF1ab_start = ()
        self.LeaderStart = 54
        self.LeaderEnd = 85
        self.AfterStartCodon = 60
        self.contig = contig
        self.reference = pysam.FastaFile(reference)
        self.refLength = self.reference.lengths[0]
        self.contig = self.reference.references[0]

    def find_start_condon(self,reference,contig,p5):

        for i in range(p5, -1, -1):
            sequence = reference.fetch(contig, i, i + 3)
            if sequence == "ATG":
                return i
        return -1

    def find_end_condon(self, reference, contig, unend_orf1b):

        for i in range(unend_orf1b, -1, -1):
            sequence = reference.fetch(contig, i, i + 3)
            if sequence == "TAA":
                return i
        return -1




    @staticmethod
    def get_read_segments_junctions_gxy(read):
        read_segments = {}
        final_read_segments = []
        final_sj_junctions = []
        sj_ref = []
        ref_start = read.reference_start
        query_start = 0
        read_conn_segment = {}
        # start
        pre_map = False
        pre_st = False
        st_refstart = 0
        st_refend = 0
        # for index, t in enumerate(read.cigartuples):
        #     st, length = t[0], t[1]
        #
        #     if st == 0 or st == 7 or st == 8:
        #         if length > 6:
        #             read_conn_segment[index] = 0
        for index, t in enumerate(read.cigartuples):
            st, length = t[0], t[1]
            if st == 3:
                pre_st=True
                # skipped region from reference
                st_refstart = ref_start-1
                st_refend = ref_start+length
                ref_start += length
                # if read_conn_segment.get(index-1) == 0 and read_conn_segment.get(index+1) == 0:
                #     sj_ref.append((ref_start - 1, ref_start + length, 'splice'))
                #     ref_start += length
            elif st == 0 or st == 7 or st == 8:
                if pre_map and pre_st and length >6:
                    sj_ref.append((st_refstart, st_refend, 'splice'))
                pre_st = False
                if length > 6:
                    pre_map = True
                else:
                    pre_map = False
                # alignment match, sequence match or sequence mismatch
                ref_start += length
                query_start += length
            elif st == 2:
                # deletion in the read
                ref_start += length  # 需要
                pre_map = False
                pre_st = False
            elif st == 1:
                # insertion
                query_start += length
                pre_map = False
                pre_st = False
            elif st == 4:
                # soft-clipping
                query_start += length
                pre_map = False
                pre_st = False
        return sj_ref




    def get_sj_reads(self, bam_file, mutList=[], paired=False) -> pd.DataFrame:
        bam = pysam.AlignmentFile(bam_file, "rb")
        count_5prime_3prime_pair = Counter()
        count = 0
        constrain_1ab = self.refLength
        splice_edge= {}
        # n_splice_edge={}
        base_ref={i : 0 for i in range(self.refLength)}
        print("The length of genome:",len(base_ref))
        for read in bam.fetch(self.contig, 0, self.refLength):

                sjs = self.get_read_segments_junctions_gxy(read) #只存在sjs,跨越边
                strand = not read.is_reverse
                if strand is True:
                    for v in read.blocks:
                        for index in range(v[0], v[1]):
                            base_ref[index] += 1  # 更新base_ref数组
                if len(sjs) > 0:
                    count += 1
                    count_5prime_3prime_pair.update(
                        [(p5, p3, strand) for p5, p3, entryType in sjs])
        # print(base_ref)
        print("the number of reads:",count)
        data = []
        for (p5, p3, strand), count in count_5prime_3prime_pair.items():
            data.append([p5, p3, strand, count])

        df = pd.DataFrame(data, columns=["5prime", "3prime", "strand", "count"])

        grouped_df = df.groupby(['5prime', '3prime']).sum()['count'].reset_index()
        grouped_df = grouped_df.sort_values('count', ascending=False)


        for _, row in grouped_df.iterrows():
            p5 = row['5prime']
            p3 = row['3prime']
            if (p5, p3) not in splice_edge.keys():
                splice_edge[(p5, p3)] = {'splice': row['count']}
            else:
                splice_edge[(p5, p3)]['splice'] += row['count']


        can_trans = {}
        for _, row in grouped_df.iterrows():
            max_count = row['count']
            break

        for _,row in grouped_df.iterrows():
            p5 = row["5prime"]
            p3 = row["3prime"]
            if  row["5prime"]>=50 and row["5prime"]<=85 :
            # if row["5prime"] >= 35 and row["5prime"] <= 85 :
               start = self.reference.fetch(self.contig, p3, self.refLength).find("ATG") + p3 #规范边的ATG开始
               un_end = start + 3 + min(self.reference.fetch(self.contig, start + 3, self.refLength).find("TAA"),
                                        self.reference.fetch(self.contig, start + 3, self.refLength).find("TAG"),
                                        self.reference.fetch(self.contig, start + 3, self.refLength).find("TGA")) #从起始密码子后一个开始找

               flag = True

               while(flag):

                   if (un_end-start)%3 == 0:
                       end = un_end+2 #un——end为终止密码子第一个的位置
                       flag = False
                   else:
                       un_end = un_end + 3 + min(
                           self.reference.fetch(self.contig, un_end + 3, self.refLength).find("TAA"),
                           self.reference.fetch(self.contig, un_end + 3, self.refLength).find("TAG"),
                           self.reference.fetch(self.contig, un_end + 3, self.refLength).find("TGA"))
                       flag=True


               is_overlap = False
               for t_start,t_end in can_trans.keys():
                   if (t_start <= start and end <= t_end) or (t_start >= start and end >= t_end):
                       is_overlap=True
                       break



               if end - start > 100 and start > 21000 and start < 29000 and not is_overlap and (
                    start, end) not in can_trans.keys():  # 转录本长度设定
                   can_trans[(start, end)] = {'splice': row['count']}
                   if start < constrain_1ab:
                       constrain_1ab = start

                   bool_index = (grouped_df["5prime"] == p5) & (grouped_df["3prime"] == p3)
                   grouped_df = grouped_df.drop(labels=grouped_df[bool_index].index, axis=0)

               elif end - start > 100 and start > 21000 and start < 29000 and not is_overlap:
                   can_trans[(start, end)]['splice'] += row['count']
                   bool_index = (grouped_df["5prime"] == p5) & (grouped_df["3prime"] == p3)
                   grouped_df = grouped_df.drop(labels=grouped_df[bool_index].index, axis=0)






        # print(can_trans)



        #print(constrain_1ab)
        constraint_1ab = constrain_1ab
        orf1ab_trans={}
        start_orf1ab = self.reference.fetch(self.contig, 200, constraint_1ab).find("ATG") + 200
        unend_orf1a= start_orf1ab + 3 + min(self.reference.fetch(self.contig, start_orf1ab + 3, constraint_1ab).find("TAA"),
                                                     self.reference.fetch(self.contig, start_orf1ab + 3, constraint_1ab).find("TAG"),
                                                     self.reference.fetch(self.contig, start_orf1ab + 3, constraint_1ab).find("TGA"))
        flag = True

        while (flag):

                    if (unend_orf1a - start_orf1ab) % 3 == 0 and unend_orf1a-start_orf1ab>13000 :
                        end_orf1a =unend_orf1a + 2
                        flag = False
                    elif unend_orf1a+3 >= constraint_1ab:
                        flag = False
                    else:
                        unend_orf1a = unend_orf1a + 3 + min(
                            self.reference.fetch(self.contig, unend_orf1a + 3, constraint_1ab).find("TAA"),
                            self.reference.fetch(self.contig, unend_orf1a + 3, constraint_1ab).find("TAG"),
                            self.reference.fetch(self.contig, unend_orf1a + 3, constraint_1ab).find("TGA"))
                        flag = True

        orf1ab_trans[(start_orf1ab, end_orf1a)] = {'orf1ab':row['count']}

        unend_orf1b = self.find_end_condon(self.reference, self.contig, constraint_1ab)

        orf1ab_trans[(start_orf1ab, unend_orf1b+2)] = {'orf1ab': row['count']}
        # print(orf1ab_trans)

        for keys ,values in can_trans.items():
            for index in range(keys[0],keys[1]):
                base_ref[index] -= values['splice']
        #print(base_ref)

        non_trans= {}
        non_trans_splice={}


        for keys in can_trans.keys():
            for _, row in grouped_df.iterrows():
                p5 = row["5prime"]
                p3 = row["3prime"]
                if p5 > keys[0] and p3 < keys[1]:
                    start = self.find_start_condon(self.reference, self.contig, p5)
                    un_end = p3 + min(self.reference.fetch(self.contig, p3, self.refLength).find("TAA"),
                                      self.reference.fetch(self.contig, p3, self.refLength).find("TAG"),
                                      self.reference.fetch(self.contig, p3, self.refLength).find("TGA"))
                    flag = True
                    while (flag):

                        if (p5 - start + un_end - p3 + 2) % 3 == 0:
                            end = un_end + 2  # un——end为终止密码子第一个的位置
                            flag = False
                        elif un_end + 3 >= 29903:
                            break
                        else:
                            un_end = un_end + 3 + min(
                                self.reference.fetch(self.contig, un_end + 3, self.refLength).find("TAA"),
                                self.reference.fetch(self.contig, un_end + 3, self.refLength).find("TAG"),
                                self.reference.fetch(self.contig, un_end + 3, self.refLength).find("TGA"))
                            flag = True
                    if tuple([(start, p5), (p3, end)]) not in non_trans.keys():
                        non_trans[(start, p5), (p3, end)] = {'non': row['count']}
                        bool_index = (grouped_df["5prime"] == p5) & (grouped_df["3prime"] == p3)
                        grouped_df = grouped_df.drop(labels=grouped_df[bool_index].index, axis=0)


                    else:
                        non_trans[(start, p5), (p3, end)]['non'] += row['count']
                        bool_index = (grouped_df["5prime"] == p5) & (grouped_df["3prime"] == p3)
                        grouped_df = grouped_df.drop(labels=grouped_df[bool_index].index, axis=0)


        for _, row in grouped_df.iterrows():
            p5 = row["5prime"]
            p3 = row["3prime"]
            if 86 < p5 < unend_orf1b:
                if p3 - p5 > 5000:
                    # start = max (self.reference.fetch(self.contig, first_values+2, p5).find("ATG") + first_values+2 )
                    # start= self.find_start_condon(self.reference,self.contig,p5)
                    start= self.find_start_condon(self.reference,self.contig,p5)
                    un_end = p3 + min(self.reference.fetch(self.contig, p3, self.refLength).find("TAA"),
                                      self.reference.fetch(self.contig, p3, self.refLength).find("TAG"),
                                      self.reference.fetch(self.contig, p3, self.refLength).find("TGA"))
                    flag = True
                    while (flag):

                        if (p5 - start + un_end - p3 + 2) % 3 == 0:
                            end = un_end + 2  # un——end为终止密码子第一个的位置
                            flag = False
                        else:
                            un_end = un_end + 3 + min(
                                self.reference.fetch(self.contig, un_end + 3, self.refLength).find("TAA"),
                                self.reference.fetch(self.contig, un_end + 3, self.refLength).find("TAG"),
                                self.reference.fetch(self.contig, un_end + 3, self.refLength).find("TGA"))
                            flag = True
                    if tuple([(start,p5),( p3,end)] )not in non_trans.keys():
                        non_trans[(start,p5),( p3,end)] = {'non': row['count']}
                        bool_index = (grouped_df["5prime"] == p5) & (grouped_df["3prime"] == p3)
                        grouped_df = grouped_df.drop(labels=grouped_df[bool_index].index, axis=0)


                    else:
                        non_trans[(start,p5),( p3,end)]['non'] += row['count']
                        bool_index = (grouped_df["5prime"] == p5) & (grouped_df["3prime"] == p3)
                        grouped_df = grouped_df.drop(labels=grouped_df[bool_index].index, axis=0)


        for t_start, t_end in can_trans.keys():
            for _, row in grouped_df.iterrows():
                p5 = row["5prime"]
                p3 = row["3prime"]
                if row["5prime"] >= 50 and row["5prime"] <= 85:
                    start = self.reference.fetch(self.contig, p3, self.refLength).find("ATG") + p3
                    un_end = start + 3 + min(
                        self.reference.fetch(self.contig, start + 3, self.refLength).find("TAA"),
                        self.reference.fetch(self.contig, start + 3, self.refLength).find("TAG"),
                        self.reference.fetch(self.contig, start + 3, self.refLength).find("TGA"))
                flag = True
                while (flag):

                    if (un_end - start) % 3 == 0:
                        end = un_end + 2  # un——end为终止密码子第一个的位置
                        flag = False
                    else:
                        un_end = un_end + 3 + min(
                            self.reference.fetch(self.contig, un_end + 3, self.refLength).find("TAA"),
                            self.reference.fetch(self.contig, un_end + 3, self.refLength).find("TAG"),
                            self.reference.fetch(self.contig, un_end + 3, self.refLength).find("TGA"))
                        flag = True

                if end > t_start and end < t_end:
                    key = tuple([(49, 84), (start, end)])
                    if key in non_trans:
                        non_trans[key]['non'] = non_trans[key].get((t_start, t_end), 0) + row['count']
                        bool_index = (grouped_df["5prime"] == p5) & (grouped_df["3prime"] == p3)
                        grouped_df = grouped_df.drop(labels=grouped_df[bool_index].index, axis=0)
                else:
                    non_trans[(49, 84), (start, end)] = {'non': row['count']}
                    bool_index = (grouped_df["5prime"] == p5) & (grouped_df["3prime"] == p3)
                    grouped_df = grouped_df.drop(labels=grouped_df[bool_index].index, axis=0)






        non_trans_after = {k: v for k, v in non_trans.items() if v.get('non', 0) >20}
        # print(non_trans_after)
        for (key1,key2),value in non_trans_after.items():
            new_key=(tuple(key1)[1],tuple(key2)[0])
            non_trans_splice[new_key]=value


        # print(non_trans_splice)



        return grouped_df, base_ref, can_trans, non_trans_after, non_trans_splice , orf1ab_trans





