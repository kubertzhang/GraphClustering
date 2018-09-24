def washImage():
    g_images = set()
    k_images = set()

    check_val = -1
    with open('rel_image', 'r') as ri_in:
        for ri_line in ri_in:
            ri_images = ri_line.strip().split(' ')
            if int(ri_images[0]) != check_val:
                check_val = int(ri_images[0])
                k_images.add(ri_images[0])
        ri_in.close()

    check_val = -1
    count = 0
    with open('rel_tag', 'r') as tag_in:
        for line in tag_in:
            img, usr, tag = line.strip().split(' ')
            if int(img) != check_val and img in k_images:
                check_val = int(img)
                g_images.add(img)
                count += 1
                if count % 20000 == 0:
                    print('loadImage = %d' % count)
        tag_in.close()

    with open('washed_img_img', 'w') as wi_ou:
        with open('rel_image', 'r') as img_in:
            for line in img_in:
                p_images = line.strip().split(' ')
                if p_images[0] not in g_images:
                    continue
                wi_ou.write('%s' % p_images[0])

                for i in range(1, len(p_images)):
                    if p_images[i] in g_images:
                        wi_ou.write('\t')
                        wi_ou.write('%s' % p_images[i])
                wi_ou.write('\n')


def imageRelationSymmetrization():
    img_img = {}
    with open('washed_img_img', 'r') as wii_in:
        for line in wii_in:
            i_imgs = line.strip().split('\t')
            img_img[i_imgs[0]] = i_imgs[1:len(i_imgs)]
    wii_in.close()

    count = 0
    for check_img in img_img:
        count += 1
        if count % 20000 == 0:
            print('syncImage = %d' % count)
        delete_list = []
        neighbor_imgs = img_img[check_img]
        for neighbor_img in neighbor_imgs:
            if check_img not in img_img[neighbor_img]:
                delete_list.append(neighbor_img)
        for de in delete_list:
            img_img[check_img].remove(de)

    with open('synced_img_img', 'w') as ii_ou:
        for key in img_img:
            ii_ou.write('%s' % key)
            for neighbor_img in img_img[key]:
                ii_ou.write('\t%s' % neighbor_img)
            ii_ou.write('\n')

# Filter data for specified data size
def washii():
    delset = set()
    with open('sii', 'r') as sii_ou:
        for line in sii_ou:
            temp_imgs = line.strip().split('\t')
            if len(temp_imgs) < 8:
                for t in temp_imgs:
                    delset.add(t)

    with open('sii2', 'w') as sii_ou:
        with open('sii', 'r') as sii2_in:
            for line in sii2_in:
                sii_imgs = line.strip().split('\t')
                if sii_imgs[0] in delset:
                    continue
                sii_ou.write('%s' % sii_imgs[0])
                for si in sii_imgs[1:len(sii_imgs)]:
                    if si in delset:
                        continue
                    sii_ou.write('\t%s' % si)
                sii_ou.write('\n')


def generateVertices():
    img_set = set()
    img_id_mapper = {}
    img_neighbor_num_mapper = {}
    img_neighbor_list_mapper = {}

    user_set = set()
    user_id_mapper = {}
    user_neighbor_num_mapper = {}
    user_neighbor_list_mapper = {}

    tag_set = set()
    tag_id_mapper = {}
    tag_neighbor_num_mapper = {}
    tag_neighbor_list_mapper = {}

    vertexid = 0

    with open('sii', 'r') as img_in:
        for line in img_in:
            tt_imgs = line.strip().split('\t')
            img_set.add(tt_imgs[0])
            img_id_mapper[tt_imgs[0]] = vertexid
            img_neighbor_num_mapper[vertexid] = [0, 0, 0]
            img_neighbor_list_mapper[vertexid] = []
            vertexid += 1

    with open('sii', 'r') as img_in2:
        for line in img_in2:
            tt_imgs2 = line.strip().split('\t')
            current_img = img_id_mapper[tt_imgs2[0]]
            img_neighbor_num_mapper[current_img][0] = len(tt_imgs2) - 1
            for t in range(1, len(tt_imgs2)):
                img_neighbor_list_mapper[current_img].append(img_id_mapper[tt_imgs2[t]])

    img_num = vertexid
    print ('img = %d' % vertexid)

    with open('rel_tag', 'r') as rtag_in:
        for line in rtag_in:
            t_img, t_user, t_tag = line.strip().split(' ')
            if t_img in img_set:
                current_img = img_id_mapper[t_img]
                if t_user in user_set:
                    current_user = user_id_mapper[t_user]
                    if current_img in user_neighbor_list_mapper[current_user]:
                        continue
                    else:
                        user_neighbor_num_mapper[current_user][0] += 1
                        user_neighbor_list_mapper[current_user].append(current_img)
                else:
                    user_set.add(t_user)
                    user_id_mapper[t_user] = vertexid
                    user_neighbor_num_mapper[vertexid] = [1, 0, 0]
                    user_neighbor_list_mapper[vertexid] = []
                    user_neighbor_list_mapper[vertexid].append(current_img)
                    vertexid += 1
                if user_id_mapper[t_user] in img_neighbor_list_mapper[current_img]:
                    continue
                else:
                    img_neighbor_num_mapper[current_img][1] += 1
                    img_neighbor_list_mapper[current_img].append(user_id_mapper[t_user])

    user_num = vertexid
    print ('user = %d' % vertexid)

    with open('rel_tag', 'r') as rtag_in2:
        for line in rtag_in2:
            t_img, t_user, t_tag = line.strip().split(' ')
            if t_img in img_set:
                current_img = img_id_mapper[t_img]
                if t_tag in tag_set:
                    current_tag = tag_id_mapper[t_tag]
                    if current_img in tag_neighbor_list_mapper[current_tag]:
                        continue
                    else:
                        tag_neighbor_num_mapper[current_tag][0] += 1
                        tag_neighbor_list_mapper[current_tag].append(current_img)
                else:
                    tag_set.add(t_tag)
                    tag_id_mapper[t_tag] = vertexid
                    tag_neighbor_num_mapper[vertexid] = [1, 0, 0]
                    tag_neighbor_list_mapper[vertexid] = []
                    tag_neighbor_list_mapper[vertexid].append(current_img)
                    vertexid += 1
                if tag_id_mapper[t_tag] in img_neighbor_list_mapper[current_img]:
                    continue
                else:
                    img_neighbor_num_mapper[current_img][2] += 1
                    img_neighbor_list_mapper[current_img].append(tag_id_mapper[t_tag])

    tag_num = vertexid
    print ('tag = %d' % vertexid)

    with open('vertices.txt', 'w') as v_ou:
        for i in range(0, img_num):
            v_ou.write('%s\t' % '\t'.join([str(i), str('image'), str(0), str(sum(img_neighbor_num_mapper[i]))]))
            v_ou.write('%s\t' % '\t'.join([str(j) for j in img_neighbor_num_mapper[i]]))
            v_ou.write('%s\n' % '\t'.join([str(k) for k in img_neighbor_list_mapper[i]]))
        for i in range(img_num, user_num):
            v_ou.write('%s\t' % '\t'.join([str(i), str('user'), str(1), str(sum(user_neighbor_num_mapper[i]))]))
            v_ou.write('%s\t' % '\t'.join([str(j) for j in user_neighbor_num_mapper[i]]))
            v_ou.write('%s\n' % '\t'.join([str(k) for k in user_neighbor_list_mapper[i]]))
        for i in range(user_num, tag_num):
            v_ou.write('%s\t' % '\t'.join([str(i), str('tag'), str(2), str(sum(tag_neighbor_num_mapper[i]))]))
            v_ou.write('%s\t' % '\t'.join([str(j) for j in tag_neighbor_num_mapper[i]]))
            v_ou.write('%s\n' % '\t'.join([str(k) for k in tag_neighbor_list_mapper[i]]))


if __name__ == '__main__':
    # washImage()
    # imageRelationSymmetrization()
    # washii()
    generateVertices()
