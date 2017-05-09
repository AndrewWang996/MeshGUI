function getAeqTest()
    meshname = 'red_dragon';
    
    Aeq = getAeq(meshname, [2, 3]);
    display(size(Aeq))
end