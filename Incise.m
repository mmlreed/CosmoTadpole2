% Static channel eroder
function g = Incise(p, g)
g.S = g.S - g.Channels.*p.E.*p.dt;

g.U(g.Channels == 1) = g.S(g.Channels == 1) - g.H(g.Channels == 1);

B = g.U;
C = g.Channels;
S = g.S;

for i=2:p.K-1
        for j=2:p.J-1
            if C(i,j)==1 & B(i,j)>=S(i+1,j) & C(i+1,j)~=1
                B(i,j)=S(i+1,j);
            elseif C(i,j)==1 & B(i,j)>=S(i,j+1) & C(i,j+1)~=1
                B(i,j)=S(i,j+1);
            elseif C(i,j)==1 & B(i,j)>=S(i,j-1) & C(i,j-1)~=1
                B(i,j)=S(i,j-1);
            elseif C(i,j)==1 & B(i,j)>=S(i-1,j-1) & C(i-1,j-1)~=1
                B(i,j)=S(i-1,j-1);
            elseif C(i,j)==1 & B(i,j)>=S(i-1,j+1) & C(i-1,j+1)~=1
                B(i,j)=S(i-1,j+1);
            elseif C(i,j)==1 & B(i,j)>=S(i+1,j+1) & C(i+1,j+1)~=1
                B(i,j)=S(i+1,j+1);
            elseif C(i,j)==1 & B(i,j)>=S(i+1,j-1) & C(i+1,j-1)~=1
                B(i,j)=S(i+1,j-1);
            elseif C(i,j)==1 & B(i,j)>=S(i-1,j) & C(i-1,j)~=1
                B(i,j)=S(i-1,j);
            end
        end
end

g.U(g.Channels==1 & g.U < p.lake_level) = p.lake_level;    

end

