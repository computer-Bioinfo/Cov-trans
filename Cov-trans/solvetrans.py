
import os
import math



class MILP_solve():

    def __init__(self, graph=None, numPaths=1,  threads=1, timelimit=None, ref_length=29903):

        self.graph = graph
        self.numPaths = numPaths
        self.threads = threads
        self.ref_length = ref_length

    def read_results(results_file):
        with open(results_file, 'r') as file:
            line = file.readline().strip()
            results = [int(float(value)) for value in line.split(',')]
        return results


    def solve_milp(self,edges,Path_edge2,Path_edge1,Path_node1,Path_node2,all_calpaths,all_nonpaths,Path,Node_edge,node1_flow,edge_flow,outfile):
        out_dir = ""
        out = open(outfile, "w")

        print('from __future__ import print_function\n', file=out, end='')
        print('import gurobipy as gp\n',file=out,end='')
        print('import sys\n', file=out, end='')
        print('import re\n', file=out, end='')
        print('def mmilp_solver(outfile):\n', file=out, end='')
        print('\ttrans = gp.Model(\"mmip\")\n',file=out,  end='')


        #Var
        # trans = gp.Model("milp_transcripts")

        for i in range(len(all_nonpaths)):
            # x=trans.addVar(lb=0,vtype=gp.GRB.BINARY,name="x_"+str(i))
            print('\tx' + str(i) + ' = trans.addVar(vtype= gp.GRB.BINARY)\n', file=out, end='')


        for i in range(len(self.graph.nodes)):
            # for j in range(i,len(self.graph.nodes)):
            #     # y = trans.addVar(lb=0,vtype=gp.GRB.BINARY, name="y_" + str(i))
            #     print('\ty' + str(i) +str(j) + ' = trans.addVar(vtype = GRB.BINARY)\n', file=out, end='')
            print('\ty' + str(i) + ' = trans.addVar(vtype = gp.GRB.BINARY)\n', file=out, end='')



        for i in range(len(edges)):
            # for j in range(i, self.graph.numEdges):
                # z = trans.addVar(lb=0,vtype=gp.GRB.BINARY, name="z_" + str(i))
                # print('\tz' + str(i) +str(j) + ' = trans.addVar(vtype=GRB.BINARY)\n', file=out, end='')
                print('\tz' + str(i) + ' = trans.addVar(vtype=gp.GRB.BINARY)\n', file=out, end='')



        for i in range(len(all_calpaths)):
            # f1 = trans.addVar(lb=0, vtype=gp.GRB.CONTINUOUS, name="f")
            print('\tf'+str(1) + str(i) + ' = trans.addVar(vtype=gp.GRB.CONTINUOUS,lb=1)\n', file=out, end='')

        for i in range(len(all_nonpaths)):
            # f2 = trans.addVar(lb=0, vtype=gp.GRB.CONTINUOUS, name="f")
            print('\tf'+str(2) + str(i) + ' = trans.addVar(vtype=gp.GRB.CONTINUOUS,lb=0)\n', file=out, end='')




        #obj
        w1=0.6
        w2=5
        w3=5

        print('\ttrans.setObjective(' + str(w1) + '*(',file=out, end='')

        for i in range(len(all_nonpaths)):
            if i < len(all_nonpaths) - 1:
                print('x' + str(i) + '+', file=out, end='')
            else:
                print('x' + str(i) + ')+' + str(w2) + '*(',file=out,end='')

        for i in range(len(self.graph.nodes)):
            if i < len(self.graph.nodes) - 1:
                print('y' + str(i) + '+', file=out, end='')
            else:
                print('y' + str(i) + ')+'+ str(w3) + '*(',file=out, end='')

        for i in range(len(edges)):
            if i < len(edges) - 1:
                print('z' + str(i) + '+', file=out, end='')
            else:
                print('z' + str(i) + '),gp.GRB.MINIMIZE)\n',file=out, end='')


        #constraint
        m=1000000
        for i in range(len(all_nonpaths)):
            # trans.addConstr(f <= m * x)
            # trans.addConstr(x <= m * f)
            print('\ttrans.addConstr(x' + str(i) + ' <= 1000000*f'+ str(2) + str(i) + ')\n',file=out, end='')
            print('\ttrans.addConstr(f' +str(2) +str(i) + ' <= 1000000*x' + str(i) + ')\n',file=out,  end='')

        # NODE constraint
        for j in range(len(self.graph.nodes)):
            temp_str = ''
            tag = 0
            temp = 0
            for i in range(len(all_calpaths)):
                temp += Path_node2[i][j]
                # print("temp_sum",temp)
            if temp > 0:
                continue
            for k in range(len(all_nonpaths)):
                temp1=Path_node1[k][j]
                if temp1 > 0:
                    if tag == 0:
                        temp_str = temp_str + 'x' + str(k)
                    else:
                        temp_str = temp_str + ' + ' + 'x' + str(k)
                    tag = tag + 1
            if temp_str != '':
                print('\ttrans.addConstr(' + temp_str + ' >= 1)\n', file=out,end='')

        # edge constraint

        for j in range(len(edges)):
            temp_str = ''
            tag = 0
            temp = 0
            for i in range(len(all_calpaths)):
                temp += Path_edge2[i][j]
            if temp > 0:
                continue
            for k in range(len(all_nonpaths)):
                temp1 = Path_edge1[k][j]
                if temp1 > 0:
                    if tag == 0:
                        temp_str = temp_str + 'x' + str(k)
                    else:
                        temp_str = temp_str + ' + ' + 'x' + str(k)
                    tag = tag + 1
            if temp_str != '':
                print('\ttrans.addConstr(' + temp_str + ' >= 1)\n', file=out,end='')

        # 点的流量总和和权重比值
        u = 0.1
        for j in range(len(self.graph.nodes)):
            temp_str1 = ''
            temp_str = ''
            tag1 = 0
            tag = 0
            for i in range(len(all_calpaths)):
                if Path_node2[i][j] == 0:
                    continue
                if tag1 == 0:
                    temp_str1 = str(Path_node2[i][j]) + '*f' + str(1)+str(i)
                    tag1 = tag1 + 1
                else:
                    temp_str1 = temp_str1 + ' + ' + str(Path_node2[i][j]) + '*f' + str(1) + str(i)

            if temp_str1 == '':
                continue

            for k in range(len(all_nonpaths)):
                if Path_node1[k][j] == 0:
                    continue
                if tag == 0:
                    temp_str = str(Path_node1[k][j]) + '*f' + str(2)+str(k)
                    tag = tag + 1
                else:
                    temp_str = temp_str + ' + ' + str(Path_node1[k][j]) + '*f' + str(2) + str(k)

            if temp_str == '':
                continue
            # print('\ttrans.addConstr(m*y' + str(j) + str(j) + ' >= ' + temp_str + '+' + str(temp) + '-' + str(
            #     (1 + u) * node1_flow[0][j]) + ')\n',file=out, end='')
            # print('\ttrans.addConstr(m*y' + str(j) + str(j) + ' >= ' + str(
            #     (1 - u) * node1_flow[0][j]) + '-'+str(temp)+' - (' + temp_str + '))\n', file=out,end='')
            print('\ttrans.addConstr(1000000*y' + str(j) + ' >= ' + temp_str + '+' + temp_str1 + '-' + str(
                (1 + 0.1) * node1_flow[0][j]) + ')\n', file=out, end='')
            print('\ttrans.addConstr(1000000*y' + str(j) + ' >= ' + str(
                (1 - 0.1) * node1_flow[0][j]) + '-(' + temp_str1 + ') - (' + temp_str + '))\n', file=out, end='')



        #边的流量总和和权重比值进行限制
        for j in range(len(edges)):
            temp_str1 = ''
            temp_str = ''
            tag1 = 0
            tag = 0
            for i in range(len(all_calpaths)):
                if Path_edge2[i][j] == 0:
                    continue
                if tag1 == 0:
                    temp_str1 = str(Path_edge2[i][j]) + '*f' + str(1) + str(i)
                    tag1 = tag1 + 1
                else:
                    temp_str1 = temp_str1 + ' + ' + str(Path_edge2[i][j]) + '*f' + str(1) + str(i)

            if temp_str1 == '':
                continue
            for k in range(len(all_nonpaths)):
                if Path_edge1[k][j] == 0:
                    continue
                if tag == 0:
                    temp_str = str(Path_edge1[k][j]) + '*f' + str(2)+str(k)
                    tag = tag + 1
                else:
                    temp_str = temp_str + ' + ' + str(Path_edge1[k][j]) + '*f' + str(2) + str(k)

            if temp_str == '':
                continue
            # print('\tmodel.addConstr(m*z' + str(j) + str(j) + ' >= ' + temp_str + '+' + str(temp) + ' - ' + str(
            #     (1 + u) * edge_flow[0][j]) + ')\n',file=out, end='')
            # print('\tmodel.addConstr(m*z' + str(j) + str(j) + ' >= ' + str(
            #     (1 - u) * edge_flow[0][j]) + '-'+str(temp) + ' - (' + temp_str + '))\n', file=out,end='')
            print('\ttrans.addConstr(1000000*z' + str(j) + ' >= ' + temp_str + '+' + temp_str1 + ' - ' + str(
                (1 + 0.1) * edge_flow[0][j]) + ')\n', file=out, end='')
            print('\ttrans.addConstr(1000000*z' + str(j) + ' >= ' + str(
                (1 - 0.1) * edge_flow[0][j]) + '- (' + temp_str1 + ') - (' + temp_str + '))\n', file=out, end='')


        print('\ttrans.optimize()\n',file=out,end='')


        print('\tofid = open(outfile, \'w\')\n', file=out, end='')

        for i in range(len(all_nonpaths)):
            if i < len(all_nonpaths) - 1:
                print('\tprint(str(x' + str(i) + '.x)+\',\', file=ofid, end=\'\')\n', file=out, end='')
            else:
                print('\tprint(str(x' + str(i) + '.x)+\'\\n\', file=ofid, end=\'\')\n', file=out, end='')
        for i in range (len(all_nonpaths)):
            if i < len(all_nonpaths)-1:
                print('\tprint(str(f2' + str(i) + '.x)+\',\', file=ofid, end=\'\')\n', file=out, end='')
            else:
                print('\tprint(str(f2' + str(i) + '.x)+\'\\n\', file=ofid, end=\'\')\n', file=out, end='')

        for i in range (len(all_calpaths)):
            if i < len(all_calpaths)-1:
                print('\tprint(str(f1' + str(i) + '.x)+\',\', file=ofid, end=\'\')\n', file=out, end='')
            else:
                print('\tprint(str(f1' + str(i) + '.x)+\'\\n\', file=ofid, end=\'\')\n', file=out, end='')

        print('mmilp_solver(\'' + out_dir + 'mmilp.results\')\n', file=out, end='')
        out.close()
        command = 'python ' + outfile + '>out_dir' + 'cov-gurobi.results'
        os.system(command)
        results_file = out_dir + "mmilp.results"
        nonresults = MILP_solve.read_results(results_file)
        return nonresults






    def writecalgtf(self, filename,transcripts_all, contig='NC_045512.2',type='cal'):
        with open(filename, 'w') as output:
            def write_transcript_and_exons(contig, path, type):
                if not path or not isinstance(path[0], (list, tuple)) or not isinstance(path[-1], (list, tuple)):
                    print(f"Invalid path: {path}")
                    return

                output.write(
                    f"{contig}\tCov-rans\ttranscript\t{path[0][0]+1}\t{path[-1][1]+1}\t.\t.\t.\t{type}\n")

                exon_start = path[0][0]
                exon_end = path[0][1]
                last_end = exon_end

                for i in range(1, len(path)):
                    if path[i][0] == last_end:
                        last_end = path[i][1]
                    else:

                        output.write(f"{contig}\tCov-rans\texon\t{exon_start+1}\t{exon_end+1}\t.\t.\t.\n")
                        exon_start = path[i][0]
                        exon_end = path[i][1]
                        last_end = exon_end


                output.write(f"{contig}\tCov-rans\texon\t{exon_start+1}\t{last_end+1}\t.\t.\t.\n")


            for transcript_info in transcripts_all:
                contig = transcript_info.get("contig", contig)
                path = transcript_info.get("path", [])
                type = transcript_info.get("type", type)
                write_transcript_and_exons(contig, path, type)

    def writenongtf(self, filename, transcripts_all, results, contig='NC_045512.2', type='non'):
        # print("Transcripts received:", transcripts_all)

        with open(filename, 'w') as output:
            def write_transcript_and_exons(contig, path):
                if not path or not isinstance(path[0], (list, tuple)) or not isinstance(path[-1], (list, tuple)):
                    print(f"Invalid path: {path}")
                    return


                output.write(f"{contig}\tCov-rans\ttranscript\t{path[0][0]+1}\t{path[-1][1]}\t.\t.\t.\t\n")


                exon_start = path[0][0]

                for i in range(1, len(path)):

                    if path[i][0] != path[i - 1][1]:

                        output.write(f"{contig}\tCov-rans\texon\t{exon_start+1}\t{path[i - 1][1]+1}\t.\t.\t.\t\n")

                        exon_start = path[i][0]



                output.write(f"{contig}\tCov-rans\texon\t{exon_start+1}\t{path[-1][1]}\t.\t.\t.\t\n")


            for transcript_info, result in zip(transcripts_all, results):
                if result == 1:
                    contig = transcript_info.get("contig", contig)
                    path = transcript_info.get("path", [])
                    write_transcript_and_exons(contig, path)
                    # print("Transcript info:", transcript_info)
                    # print("Path:", path)
        
