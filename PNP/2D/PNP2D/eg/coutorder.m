clear
clc


NO=0;  % Please change this number for the following case.


switch NO

    case 0
Error=[ 
1.0673e-06 2.2235e-07 4.9358e-07
1.3252e-07 2.7599e-08 5.7668e-08
3.9247e-08 8.1674e-09 1.6726e-08
1.6791e-08 3.502e-09 7.0031e-09
]';

    case 1
Error=[ 
0.13404  0.0451084
0.0381464  0.0124465
0.0100183  0.0031893
0.00255603  0.000802206

]';
    case 2
Error=[
0.00740928  0.00322366
0.000991089  0.000435251
0.000126183  5.55145e-05
1.58469e-05  6.97483e-06
]';        
    case 3
Error=[
0.000553341  0.000360523
3.77612e-05  2.38081e-05
2.41654e-06  1.51433e-06
1.51952e-07  9.51458e-08
]';
    case 4
Error=[ 

]';       
    otherwise
        disp('a wrong input!')
end

[m,n]=size(Error);
order=zeros(m,n-1);
for i=1:m
    for j=1:n-1
        order(i,j)=log(Error(i,j)/Error(i,j+1))/log((j+1)/j);
    end
end

order=order';
[M,N]=size(order);
for i=1:M
    fprintf('%7.2f',order(i,:))
    fprintf('\n')
end

ErrSav=Error';

%stop

if NO==0
    stream = fopen('1dp0.txt', 'w');
elseif NO==1
    stream = fopen('1dp1.txt', 'w');
elseif NO==2
    stream = fopen('1dp2.txt', 'w');
elseif NO==3
    stream = fopen('1dp3.txt', 'w');
elseif NO==4
    stream = fopen('1dp4.txt', 'w');
end

fprintf(stream, '%g\n', n);
for i=1:n
    for j=1:m
        fprintf(stream, '%g\t', ErrSav(i,j));
    end
    fprintf(stream, '\n');
end
fprintf(stream, '\n');
fprintf(stream, '%g\n', M);
for i=1:M
    fprintf(stream, '%7.2f',order(i,:));
    fprintf(stream, '\n');
end
fclose(stream);
