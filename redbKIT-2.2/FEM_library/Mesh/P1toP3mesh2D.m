function [elements, vertices, boundaries, rings] = P1toP3mesh2D(elements,vertices,boundaries, rings)
%P1TOP3MESH2D builds a P3 mesh in 2D.
%       F. Saleri 9-20-01.
%       F. Negri 2016, Add mesh Graph to speedup computations (still inefficient)
%       S. Cominelli 2022, Add P3

[~,nov]     =  size(vertices);
[~,noe]     =  size(elements);
nside       =  nov;

elements = [elements(1:3,:); zeros(7,noe); elements(4,:)];

[ a ] = compute_adjacency(vertices, elements, 2, 'P1');

[ii,jj,vv] = find(a);
a = sparse(ii, jj, vv*0 - 1, nov, nov);
pesi = [2  1; 1  2]/3;
for ie = 1:noe
    i = elements(1,ie);
    j = elements(2,ie);
    k = elements(3,ie);
    
    l1 = a(i,j);
    if l1 == -1 
        nside = nside + 1;
        a(i,j) = nside;
        a(j,i) = nside;
        elements(4,ie) = nside;
        nside = nside + 1;
        elements(5,ie) = nside;
        vertices(1,nside-1:nside) = pesi * vertices(1,[i j]).';
        vertices(2,nside-1:nside) = pesi * vertices(2,[i j]).';
    else % devo assicurare che l'odine dei nodi all'interno della colonna elements(ie,:) sia i -> j e non j -> i
        if abs(diff(vertices(1,[i j]))) > 10000*eps % se le x non sono troppo uguali, guardo la condizione sulle x
            verso_giusto = abs(vertices(1,l1) - pesi(1,:)*vertices(1,[i j]).') < 100*eps;
        else    % altrimenti guardo la condizione sulle y
            verso_giusto = abs(vertices(2,l1) - pesi(1,:)*vertices(2,[i j]).') < 100*eps;
        end
        if verso_giusto   % in sostanza if vertices(1,l1) == pesi(1,:)*vertices(1,[i j]).'
            elements(4,ie) = l1;
            elements(5,ie) = l1 + 1;
        else
            elements(4,ie) = l1 + 1;
            elements(5,ie) = l1;
        end
    end
   
    l2 = a(j,k);
    if l2 == -1 
        nside = nside + 1;
        a(j,k) = nside;
        a(k,j) = nside;
        elements(6,ie) = nside;
        nside = nside + 1;
        elements(7,ie) = nside;
        vertices(1,nside-1:nside) = pesi * vertices(1,[j k]).';
        vertices(2,nside-1:nside) = pesi * vertices(2,[j k]).';
    else 
        if abs(diff(vertices(1,[j k]))) > 10000*eps % se le x non sono troppo uguali, guardo la condizione sulle x
            verso_giusto = abs(vertices(1,l2) - pesi(1,:)*vertices(1,[j k]).') < 100*eps;
        else    % altrimenti guardo la condizione sulle y
            verso_giusto = abs(vertices(2,l2) - pesi(1,:)*vertices(2,[j k]).') < 100*eps;
        end
        if verso_giusto
            elements(6,ie) = l2;
            elements(7,ie) = l2 + 1;
        else
            elements(6,ie) = l2 + 1;
            elements(7,ie) = l2;
        end
    end
   
    l3 = a(k,i);
    if l3 == -1 
        nside = nside + 1;
        a(k,i) = nside;
        a(i,k) = nside;
        elements(8,ie) = nside;
        nside = nside + 1;
        elements(9,ie) = nside;
        vertices(1,nside-1:nside) = pesi * vertices(1,[k i]).';
        vertices(2,nside-1:nside) = pesi * vertices(2,[k i]).';
    else
        if abs(diff(vertices(1,[k i]))) > 10000*eps % se le x non sono troppo uguali, guardo la condizione sulle x
            verso_giusto = abs(vertices(1,l3) - pesi(1,:)*vertices(1,[k i]).') < 100*eps;
        else    % altrimenti guardo la condizione sulle y
            verso_giusto = abs(vertices(2,l3) - pesi(1,:)*vertices(2,[k i]).') < 100*eps;
        end
        if verso_giusto
            elements(8,ie) = l3;
            elements(9,ie) = l3 + 1;
        else
            elements(8,ie) = l3 + 1;
            elements(9,ie) = l3;
        end
    end
    
    nside = nside + 1;
    elements(10,ie) = nside;
    vertices(1,nside) = sum(vertices(1,[i j k]))/3;
    vertices(2,nside) = sum(vertices(2,[i j k]))/3;
end

[~,nside] = size(boundaries);
for i = 1 : nside
    j = boundaries(1,i);
    k = boundaries(2,i);
    l1 = a(j,k);

    if abs(diff(vertices(1,[j k]))) > 10000*eps % se le x non sono troppo uguali, guardo la condizione sulle x
        verso_giusto = abs(vertices(1,l1) - pesi(1,:)*vertices(1,[j k]).') < 100*eps;
    else    % altrimenti guardo la condizione sulle y
        verso_giusto = abs(vertices(2,l1) - pesi(1,:)*vertices(2,[j k]).') < 100*eps;
    end
    if verso_giusto   % in sostanza if vertices(1,l1) == pesi(1,:)*vertices(1,[k i]).'
        boundaries(3,i) = l1;
        boundaries(4,i) = l1 + 1;
    else
        boundaries(3,i) = l1 + 1;
        boundaries(4,i) = l1;
    end
end

end