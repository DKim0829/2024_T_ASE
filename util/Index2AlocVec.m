function aloc_vec = Index2AlocVec(index, num_total_goals)
    aloc_vec = (dec2bin(index-1, num_total_goals)=='1');
end