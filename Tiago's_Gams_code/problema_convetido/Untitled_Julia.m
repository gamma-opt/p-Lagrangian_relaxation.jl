clear all
clc
fileID = fopen('C:\\Artigos\\p-LD\\numerical experiments\\v7\\problema_convetido\\problema.csv','r');
pro = fscanf(fileID,'%s');
problema = ['C:\\Artigos\\p-LD\\numerical experiments\\problems qcqp\\qcqp\\qcqp\\' pro '.m'];
run(problema)

n = size(x0);
n = n(1);
m = size(cl);
m = m(1);

x(1:n) = 0;
a_(1:n,1:m) = 0;
a(1:n,1:n,1:m) = 0;
b(1:n,1:m) = 0;
a__FO(1:n) = 0;
a_FO(1:n,1:n) = 0;
b_FO(1:n) = 0;
c = nlcon(x);
c_FO = fun(x);

% for i = 1:n
%     for j = 1:n
%         x(j) = 0;
%     end
%     x(i) = 1;
%     a_(i,:) = nlcon(x);
% end
% 
% for i = 1:n
%     for j = 1:n
%         x(j) = 0;
%     end
%     x(i) = 1;
%     a_(i,:) = nlcon(x);
%     a__FO(i) = fun(x);
% end

for i = 1:n
    for j = 1:n
        x(j) = 0;
    end
    x(i) = 1;
    f1 = nlcon(x);
    f1_FO = fun(x);
    x(i) = -1;
    f2 = nlcon(x);
    f2_FO = fun(x);
    a_(i,:) = (f1+f2)/2 - c;
    a__FO(i) = (f1_FO+f2_FO)/2 - c_FO;
    b(i,:) = (f1-f2)/2 ;
    b_FO(i) = (f1_FO-f2_FO)/2 ;
end


aux(1,1:n) = 0;
for i = 1:n
    for j = 1:n
        x(j) = 0;
    end
    x(i) = 1;
    for k = i:n
        aux(k) = 1;
        x = x + aux;
        if i == k
            a(i,i,:) = a_(i,:);
            a_FO(i,i) = a__FO(i);
        else
            a(i,k,:) =  nlcon(x).' - b(i,:) - b(k,:) - c.' - a_(i,:) - a_(k,:);
            a_FO(i,k) =  fun(x) - b_FO(i) - b_FO(k) - c_FO - a__FO(i) - a__FO(k);
        end
        x = x - aux;
        aux(k) = 0;
    end
end


% a_ = a_.';
%b = b.';
% 
% a__FO = a__FO.';
% b_FO = b_FO.';


% Descobre autovalores da matriz
% Primeiro transforma a matriz triangular em uma matriz simétrica e depois
% pega o seu autovalor


%método 1: identidade
%autovalores da função objetivo
a_sym=a_FO'+a_FO;
a_sym(1:n+1:end) = diag(a_FO);
eigen = eig(a_sym);
eigen_min_FO = min(eigen);

%autovalores das restrições
eigen_min(1:k) = 0;
for k = 1:m
    a_sym=a(:,:,k)'+a(:,:,k);
    a_sym(1:n+1:end) = diag(a(:,:,k));
    eigen = eig(a_sym);
    eigen_min(k) = min(eigen);   
end

%método 2: diagonal dominant
a_sym=a_FO'+a_FO;
a_sym(1:n+1:end) = diag(a_FO);

r_DD_FO(i:n) = 0;
for i = 1:n
    r_DD_FO(i) = -a_sym(i,i);
    for j = 1:n
        if not(j == i)
            r_DD_FO(i) = r_DD_FO(i) + abs(a_sym(i,j));
        end
    end
    r_DD_FO(i) = max(0, r_DD_FO(i));
end

r_DD(i:n, 1:m) = 0;
for k = 1:m
    a_sym=a(:,:,k)'+a(:,:,k);
    a_sym(1:n+1:end) = diag(a(:,:,k));

    for i = 1:n
        r_DD(i,k) = -a_sym(i,i);
        for j = 1:n
            if not(j == i)
                r_DD(i,k) = r_DD(i,k) + abs(a_sym(i,j));
            end
        end
        r_DD(i,k) = max(0, r_DD(i,k));
    end
end


% método 3: Diagonal SDP
% parte da função objetivo
a_sym=a_FO'+a_FO;
a_sym(1:n+1:end) = diag(a_FO);

% cvx_begin
%     variable x(n)
%     variable P(n,n) symmetric
%     dual variable pi
%     minimize(sum(x) )
%     subject to
%         pi: x >= 0;
%         a_sym + diag(x) == semidefinite(n);
% cvx_end
% r_DSDP_FO = x;
% 
% % parte das restrições
% r_DSDP(1:n, 1:m) = 0;
% for k = 1:m
%     a_sym=a(:,:,k)'+a(:,:,k);
%     a_sym(1:n+1:end) = diag(a(:,:,k));
% 
%     cvx_begin
%         variable x(n)
%         variable P(n,n) symmetric
%         dual variable pi
%         minimize(sum(x) )
%         subject to
%             pi: x >= 0;
%             a_sym + diag(x) == semidefinite(n);
%     cvx_end
%     r_DSDP(:,k) = x;
% end


%%%%%%%%%%%%%%%%%%
% escreve em arquivo
fileID_aux = fopen('C:\\Users\\tiago.andrade\\Desktop\\RNMDT_Julia\\Generated instances\\instance.jl','w');
newText_aux{1} = ['n =' int2str(n)];
newText_aux{2} = ['m =' int2str(m)];
newText_aux{3} = ['lb = zeros(n)'];
newText_aux{4} = ['ub = ones(n)'];
newText_aux{5} = ['c = zeros(m+1)'];
newText_aux{6} = ['b = zeros(m+1,n)'];
newText_aux{7} = ['Q = zeros(m+1,n,n)'];
fprintf(fileID_aux,'%s\n',newText_aux{:});
fclose(fileID_aux);


fileID = fopen('C:\\Users\tiago.andrade\\Desktop\RNMDT_Julia\\Generated instances\\instance.jl','a');
linha = 0;
newText{1} = '';
for i=1:n
    if (b_FO(i) <= -10^-6 || b_FO(i) >= 10^-6)
        linha = linha + 1;
        newText{linha} = ['b[' '1' ',' int2str(i) '] = ' num2str(b_FO(i))]; 
    end
end

for i=1:n
    for j = 1:m
        if (b(i,j) <= -10^-6 || b(i,j) >= 10^-6)
            linha = linha + 1;
            newText{linha} = ['b[' int2str(j+1) ',' int2str(i) '] = ' num2str(b(i,j))]; 
        end
    end
end

linha = linha + 1;
newText{linha} = ['c[' '1' '] = ' num2str(-c_FO)]; 

for j = 1:m
        linha = linha + 1;
        newText{linha} = ['c[' int2str(j+1) '] = ' num2str(c(j)-cu(j))]; 
end


for i = 1:n
        linha = linha + 1;
        newText{linha} = ['lb[' int2str(i) '] = ' num2str(lb(i))]; 
        linha = linha + 1;
        newText{linha} = ['ub[' int2str(i) '] = ' num2str(ub(i))]; 
end

for i=1:n
    for j = 1:n
        for k = 1:m
            if (a_FO(i,j) <= -10^-6 || a_FO(i,j) >= 10^-6)
                %linha = linha + 1;
                %newText{linha} = ['BL(' '1' ',' int2str(i) ',' int2str(j) ') = ' '1']; 
                linha = linha + 1;
                newText{linha} = ['Q[' '1' ',' int2str(i) ',' int2str(j) '] = ' num2str(a_FO(i,j))]; 
            end
        end
    end
end
for i=1:n
    for j = 1:n
        for k = 1:m
            if (a(i,j,k) <= -10^-6 || a(i,j,k) >= 10^-6)
                %linha = linha + 1;
                %newText{linha} = ['BL(' int2str(k+1) ',' int2str(i) ',' int2str(j) ') = ' '1']; 
                linha = linha + 1;
                newText{linha} = ['Q[' int2str(k+1) ',' int2str(i) ',' int2str(j) '] = ' num2str(a(i,j,k))]; 
            end
        end
    end
end

fprintf(fileID,'%s\n',newText{:});
fclose(fileID);


% fileID = fopen('C:\\Artigos\\p-LD\\numerical experiments\\v6\\problema_convetido\\r_DD.inc','w');
% linha = 0;
% newText_r_dd{1} = '';
% linha = 0;
% newText_r_dd{1} = '';
% for i = 1:n
%     linha = linha + 1;
%     newText_r_dd{linha} = ['r_DD(' char(39) '0' char(39) ',' char(39) int2str(i) char(39) ') = ' num2str(r_DD_FO(i))]; 
% end
% for k = 1:m
%     for i=1:n
%             linha = linha + 1;
%             newText_r_dd{linha} = ['r_DD(' char(39) int2str(k) char(39) ',' char(39) int2str(i) char(39) ') = ' num2str(r_DD(i,k))]; 
%     end
% end
% fprintf(fileID,'%s\n',newText_r_dd{:});
% fclose(fileID);
% 
% 
% fileID = fopen('C:\\Artigos\\p-LD\\numerical experiments\\v6\\problema_convetido\\r_DSDP.inc','w');
% linha = 0;
% newText_r_dsdp{1} = '';
% linha = 0;
% newText_r_dsdp{1} = '';
% for i = 1:n
%     linha = linha + 1;
%     newText_r_dsdp{linha} = ['r_dsdp(' char(39) '0' char(39) ',' char(39) int2str(i) char(39) ') = ' num2str(r_DSDP_FO(i))]; 
% end
% for k = 1:m
%     for i=1:n
%             linha = linha + 1;
%             newText_r_dsdp{linha} = ['r_dsdp(' char(39) int2str(k) char(39) ',' char(39) int2str(i) char(39) ') = ' num2str(r_DSDP(i,k))]; 
%     end
% end
% fprintf(fileID,'%s\n',newText_r_dsdp{:});
% fclose(fileID);