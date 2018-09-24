
with open('football_vertices', 'r') as dbou:
    checklist = []
    for line in dbou:
        checklist = line.strip().split('\t')
        for id in checklist[7:len(checklist)]:
            if id == checklist[0]:
                print(id)
        for i in range(7, len(checklist)):
            if checklist[i] in checklist[i+1 : len(checklist)]:
                    print(checklist[0])