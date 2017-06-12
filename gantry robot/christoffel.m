function C = christoffel(m_matrix)
syms tht2 tht3 tht4 L1 L2 L3 L4
vars = [L1 tht2 tht3 tht4]

for i = 1:4
    for j=1:4
        
        for k = 0:4
            bijk = 1/2*(diff(m_matrix(i,j),vars(k)) + diff(m_matrix(i,k),vars(j)) - diff(m_matrix(j,k),vars(i)))
        end      
        
    end
end

end
