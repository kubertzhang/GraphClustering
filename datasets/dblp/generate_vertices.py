# -*- coding: utf-8 -*-

conf_and_jour = set(['432E313B', '432E31B6', '432E3307', '432E2FC4'])

conf_jour_count = {}


# 读取conf_and_jour集合指定的 conference 和 journal
def loadConferenceAndJournal():
    res = set()
    for element in conf_and_jour:
        res.add(int(element, base=16))
        conf_jour_count[int(element, base=16)] = 0   # 初始化会议论文计数
    return res


# 筛选指定conf_and_jour内的paper
def wash_paper():
    no_place = 0
    with open('CleanPapers.txt', 'r') as f:
        dest = open('SmallCleanPapers.txt', 'w')
        line_num = 0
        for line in f:
            line_num += 1
            segments = line.strip().split('\t')
            paper_id = segments[0]
            place_id = segments[1]
            if place_id in conf_and_jour:
                # ---- 二选一
                #dest.write("%s\t%s\n" % (paper_id, place_id))
                #  --- 规定选择指定数目的 paper 数量 ---
                if place_id in conf_jour_count:
                    conf_jour_count[place_id] += 1
                else:
                    conf_jour_count[place_id] = 0
                    
                if conf_jour_count[place_id] < 125:        # 单个会议的paper数目达到1000篇时，不再添加该论文的paper
                    dest.write("%s\t%s\n" % (paper_id, place_id))
                # --------------------------------------
            else:
                no_place += 1
            if line_num % 200000 == 0:
                print line_num, no_place

                
# 将指定conf_and_jour内的paper加载到集合: 过滤重复的paper
def loadPapers():
    paper_input = open('NewPapers.txt', 'r')
    papers = set()
    line_cnt = 0
    for line in paper_input:
        line_cnt += 1
        papers.add(int(line.strip(), base=16))
        if line_cnt % 200000 == 0:
            print line_cnt
    return papers


def is_hex(num):
    if len(num) != 8:
        return False
    return all(['A' <= ch <= 'F' or '0' <= ch <= '9' for ch in num])


#def wash_keyword():
#    with open('FieldsOfStudy.txt', 'r') as f:
#        dest = open('CleanKeywords.txt', 'w')
#        line_num = 0
#        for line in f:
#            line_num += 1
#            segments = line.split('\t')
#            dest.write("%s\n" % (segments[0]))
#            if line_num % 100000 == 0:
#                print line_num


def wash_edges():

    p = loadPapers()

    # 处理paper - author 
    with open('CleanPaperAuthorAffiliations.txt', 'r') as painput:
        paoutput = open('SmallCleanPaperAuthorAffiliations.txt', 'w')
        
        legal_edges = 0
        line_cnt = 0
        for line in painput:
            line_cnt += 1
            paperid, authorid = line.strip().split('\t')[0:2]
            if int(paperid, base=16) in p:
                legal_edges += 1
                paoutput.write('%s\t%s\n' % (paperid, authorid))
            if line_cnt % 200000 == 0:
                print 'line_no, legal: ', line_cnt, legal_edges
        paoutput.close()

    # 处理paper - paper
    with open('CleanPaperReferences.txt', 'r') as ppinput:
        ppoutput = open('SmallCleanPaperReferences.txt', 'w')
        legal_edges = 0
        line_cnt = 0
        for line in ppinput:
            line_cnt += 1
            paperid, refid = line.strip().split('\t')[0:2]
            if int(paperid, base=16) in p and int(refid, base=16) in p and int(paperid, base=16) != int(refid, base=16):
                legal_edges += 1
                ppoutput.write('%s\t%s\n' % (paperid, refid))
                ppoutput.flush()
            if line_cnt % 200000 == 0:
                print 'line_no, legal: ', line_cnt, legal_edges
        ppoutput.close()

    # 处理paper - keywords
    with open('CleanPaperKeywords.txt', 'r') as pkinput:
        pkoutput = open('SmallCleanPaperKeywords.txt', 'w')
        legal_edges = 0
        line_cnt = 0
        for line in pkinput:
            line_cnt += 1
            paperid, keywordid = line.strip().split('\t')[0:2]
            if int(paperid, base=16) in p:
                legal_edges += 1
                pkoutput.write('%s\t%s\n' % (paperid, keywordid))
                pkoutput.flush()
            if line_cnt % 200000 == 0:
                print 'line_no, legal: ', line_cnt, legal_edges
        pkoutput.close()


def to_hex(num):
    res = ''
    for i in xrange(8):
        last = num % 16
        ch = str(last) if 0 <= last <= 9 else chr(ord('A')+last-10)
        res = ch + res
        num /= 16
    return res


def generate_vertices():

    vertex = 0

    p_id_mapper = {}
    a_id_mapper = {}
    c_id_mapper = {}
    k_id_mapper = {}

    p_edge_mapper = {}
    a_edge_mapper = {}
    c_edge_mapper = {}
    k_edge_mapper = {}

    p_neighbor_mapper = {}
    a_neighbor_mapper = {}
    c_neighbor_mapper = {}
    k_neighbor_mapper = {}

    # 过滤处理
    pp_check_table = {}
    pa_check_table = {}
    pc_check_table = {}
    pk_check_table = {}

    p = loadPapers()
    c = loadConferenceAndJournal()

    # 对paper进行编号
    for paper in p: 
        p_id_mapper[paper] = vertex                # paper 编号
        p_edge_mapper[paper] = [0, 0, 0, 0]        # paper 对应的4类邻居点数目
        p_neighbor_mapper[paper] = []              # 用来存储邻居点
        
        # 初始化筛选集合
        pp_check_table[paper] = []
        pa_check_table[paper] = []
        pc_check_table[paper] = []
        pk_check_table[paper] = []

        vertex += 1

    print 'caculating p-p relation'
    with open('SmallCleanPaperReferences.txt', 'r') as ppinput:
        for line in ppinput:
            paperid, refid = line.strip().split('\t')[0:2]
            paperid = int(paperid, base=16)
            refid = int(refid, base=16)

            ##
            if refid in pp_check_table[paperid] or paperid in pp_check_table[refid]:
                continue

            pp_check_table[paperid].append(refid)
            pp_check_table[refid].append(paperid)
            ##

            p_edge_mapper[paperid][0] += 1
            p_edge_mapper[refid][0] += 1

            p_neighbor_mapper[paperid].append(p_id_mapper[refid])
            p_neighbor_mapper[refid].append(p_id_mapper[paperid])

    print 'caculating p-a relation'
    with open('SmallCleanPaperAuthorAffiliations.txt', 'r') as painput:
        for line in painput:
            paperid, authorid = line.strip().split('\t')[0:2]
            paperid = int(paperid, base=16)
            authorid = int(authorid, base=16)

            ##
            if authorid in pa_check_table[paperid]:
                continue

            pa_check_table[paperid].append(authorid)
            ##
            
            if authorid not in a_id_mapper:
                a_id_mapper[authorid] = vertex
                a_edge_mapper[authorid] = 0
                a_neighbor_mapper[authorid] = []
                vertex += 1
            a_edge_mapper[authorid] += 1
            p_edge_mapper[paperid][1] += 1
            p_neighbor_mapper[paperid].append(a_id_mapper[authorid])
            a_neighbor_mapper[authorid].append(p_id_mapper[paperid])

    for conf in c:
        c_id_mapper[conf] = vertex
        c_edge_mapper[conf] = 0
        c_neighbor_mapper[conf] = []
        vertex += 1

    print 'caculating p-c relation'
    with open('SmallCleanPapers.txt', 'r') as pinput:
        for line in pinput:
            paperid, confid = line.strip().split('\t')[0:2]
            paperid = int(paperid, base=16)
            confid = int(confid, base=16)

            ##
            if paperid not in p or confid in pc_check_table[paperid]:
                continue

            pc_check_table[paperid].append(confid)
            ##
            
            c_edge_mapper[confid] += 1
            p_edge_mapper[paperid][2] += 1
            p_neighbor_mapper[paperid].append(c_id_mapper[confid])
            c_neighbor_mapper[confid].append(p_id_mapper[paperid])

    print 'caculating p-k relation'
    with open('SmallCleanPaperKeywords.txt', 'r') as pkinput:
        for line in pkinput:
            paperid, keywordid = line.strip().split('\t')[0:2]
            paperid = int(paperid, base=16)
            keywordid = int(keywordid, base=16)

            ##
            if keywordid in pk_check_table[paperid]:
                continue

            pk_check_table[paperid].append(keywordid)
            ##

            if keywordid not in k_id_mapper:
                k_id_mapper[keywordid] = vertex
                k_edge_mapper[keywordid] = 0
                k_neighbor_mapper[keywordid] = []
                vertex += 1
            k_edge_mapper[keywordid] += 1
            p_edge_mapper[paperid][3] += 1
            p_neighbor_mapper[paperid].append(k_id_mapper[keywordid])
            k_neighbor_mapper[keywordid].append(p_id_mapper[paperid])

    with open('SmallVertices.txt', 'w') as vout:
        for paper in p:
            vout.write('%s\t%s\t0\t%s\t%s\t%s\n' % (p_id_mapper[paper], to_hex(paper), str(sum(p_edge_mapper[paper])),
                                                    '\t'.join([str(i) for i in p_edge_mapper[paper]]),
                                                    '\t'.join([str(j) for j in p_neighbor_mapper[paper]])))

        authors = sorted(a_id_mapper.items(), lambda x, y: cmp(x[1], y[1]))
        for author in authors:
            vout.write('%s\t%s\t1\t%s\t%s\t%s\n' % (author[1], to_hex(author[0]), a_edge_mapper[author[0]],
                                                    '\t'.join([str(a_edge_mapper[author[0]]), '0', '0', '0']),
                                                    '\t'.join([str(j) for j in a_neighbor_mapper[author[0]]])))

        for conf in c:
            vout.write('%s\t%s\t2\t%s\t%s\t%s\n' % (c_id_mapper[conf], to_hex(conf), c_edge_mapper[conf],
                                                    '\t'.join([str(c_edge_mapper[conf]), '0', '0', '0']),
                                                    '\t'.join([str(j) for j in c_neighbor_mapper[conf]])))

        keywords = sorted(k_id_mapper.items(), lambda x, y: cmp(x[1], y[1]))
        for key in keywords:
            vout.write('%s\t%s\t3\t%s\t%s\t%s\n' % (key[1], to_hex(key[0]), k_edge_mapper[key[0]],
                                                    '\t'.join([str(k_edge_mapper[key[0]]), '0', '0', '0']),
                                                    '\t'.join([str(j) for j in k_neighbor_mapper[key[0]]])))


if __name__ == '__main__':
    # 筛选指定conf_and_jour, 生成小数据集
    # wash_paper()
    wash_edges()
    # 生成结果节点文件
    generate_vertices()

