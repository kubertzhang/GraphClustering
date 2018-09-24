total_neighbor_index = 3
player_neighbor_index = 4
list_neighbor_index = 5
tweet_neighbor_index = 6

players_mapper = {}
vertices_mapper = {}

cnt = 0
# load players
with open('football/football.communities', 'r') as fc_in:
    for line in fc_in:
        seg = line.strip().split(':')[1].strip().split(',')
        for player in seg:
            players_mapper[int(player)] = cnt
            vertices_mapper[cnt] = []
            vertices_mapper[cnt].extend([cnt, int(player), 0, 0, 0, 0, 0])
            cnt += 1

# load player-player relationship
skip_mark = True
with open('football/football-followedby.mtx', 'r') as ff_in:
    for line in ff_in:
        if skip_mark:
            skip_mark = False
            continue
        seg = line.strip().split(' ')
        player1 = players_mapper[int(seg[0])]
        player2 = players_mapper[int(seg[1])]

        if player2 not in vertices_mapper[player1][7:]:
            vertices_mapper[player1][total_neighbor_index] += 1
            vertices_mapper[player1][player_neighbor_index] += 1
            vertices_mapper[player1].append(player2)

        if player1 not in vertices_mapper[player2][7:]:
            vertices_mapper[player2][total_neighbor_index] += 1
            vertices_mapper[player2][player_neighbor_index] += 1
            vertices_mapper[player2].append(player1)

# load list-player relationship
skip_mark = True
tmp_cnt = cnt
with open('football/football-listmerged500.mtx', 'r') as fl_in:
    for line in fl_in:
        if skip_mark:
            skip_mark = False
            continue
        seg = line.strip().split(' ')
        list_id = tmp_cnt + int(seg[0])
        player_id = players_mapper[int(seg[1])]
        if list_id not in vertices_mapper.keys():
            cnt += 1
            vertices_mapper[list_id] = []
            vertices_mapper[list_id].extend([list_id, list_id, 1, 0, 0, 0, 0])

        vertices_mapper[list_id][total_neighbor_index] += 1
        vertices_mapper[list_id][player_neighbor_index] += 1
        vertices_mapper[list_id].append(player_id)

        vertices_mapper[player_id][total_neighbor_index] += 1
        vertices_mapper[player_id][list_neighbor_index] += 1
        vertices_mapper[player_id].append(list_id)

# load tweet-player relationship
skip_mark = True
with open('football/football-tweets500.mtx', 'r') as ft_in:
    for line in ft_in:
        if skip_mark:
            skip_mark = False
            continue
        seg = line.strip().split(' ')
        tweet_id = cnt + int(seg[0])
        player_id = players_mapper[int(seg[1])]
        if tweet_id not in vertices_mapper.keys():
            vertices_mapper[tweet_id] = []
            vertices_mapper[tweet_id].extend([tweet_id, tweet_id, 2, 0, 0, 0, 0])

        vertices_mapper[tweet_id][total_neighbor_index] += 1
        vertices_mapper[tweet_id][player_neighbor_index] += 1
        vertices_mapper[tweet_id].append(player_id)

        vertices_mapper[player_id][total_neighbor_index] += 1
        vertices_mapper[player_id][tweet_neighbor_index] += 1
        vertices_mapper[player_id].append(tweet_id)

# store result
fts = sorted(vertices_mapper.items(), cmp=lambda x, y: cmp(x[0], y[0]))
with open('football_vertices', 'w') as fv_ou:
    for ft in fts:
        fv_ou.write('%s\n' % '\t'.join([str(c) for c in ft[1]]))




