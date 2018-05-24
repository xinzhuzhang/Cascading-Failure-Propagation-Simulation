clc;
clear;

load Basedata.mat;

N1index=find(DeD==1);
N2index=find(DeD==2);
N3index=find(DeD==3); 
N4index=find(DeD==4);
N5index=find(DeD==5);
N6index=find(DeD==6);

LN1=zeros(100,10);
LN2=zeros(100,10);
LN3=zeros(100,10);
LN4=zeros(100,10);
LN5=zeros(100,10);
LN6=zeros(100,10);

NL1=zeros(100,11);
NL2=zeros(100,11);
NL3=zeros(100,11);
NL4=zeros(100,11);
NL5=zeros(100,11);
NL6=zeros(100,11);

for y=1:100 
    rescueNode4=N4index(ceil(rand()*size(N4index,2)));
    rescue_pathlength=zeros(num,size(rescueNode4,2));
    rescue_path=cell(num,size(rescueNode4,2));
    for i=1:num
        for j=1:size(rescueNode4,2)
            [rescue_pathlength(i,j) rescue_path{i,j}]=graphshortestpath(sparse(matrix_length),rescueNode4(j),node(i));
        end
    end
    rescue_meanpathlength=mean(rescue_pathlength);
    
N1index_Lessthan=[];
N2index_Lessthan=[];
N3index_Lessthan=[]; 
N4index_Lessthan=[];
N5index_Lessthan=[];
N6index_Lessthan=[];
N1index_Greaterthan=[];
N2index_Greaterthan=[];
N3index_Greaterthan=[]; 
N4index_Greaterthan=[];
N5index_Greaterthan=[];
N6index_Greaterthan=[];
Tmax=zeros(1,num);
for i=1:num
    if rescue_pathlength(i)<rescue_meanpathlength
        Tmax(i)=Tmax(i)+1;
        switch DeD(i)
            case 1
                N1index_Lessthan=[N1index_Lessthan i];
            case 2
                N2index_Lessthan=[N2index_Lessthan i];
            case 3
                N3index_Lessthan=[N3index_Lessthan i];
            case 4
                N4index_Lessthan=[N4index_Lessthan i];
            case 5
                N5index_Lessthan=[N5index_Lessthan i];
            case 6
                N6index_Lessthan=[N6index_Lessthan i];
        end
    else
        switch DeD(i)
            case 1
                N1index_Greaterthan=[N1index_Greaterthan i];
            case 2
                N2index_Greaterthan=[N2index_Greaterthan i];
            case 3
                N3index_Greaterthan=[N3index_Greaterthan i];
            case 4
                N4index_Greaterthan=[N4index_Greaterthan i];
            case 5
                N5index_Greaterthan=[N5index_Greaterthan i];
            case 6
                N6index_Greaterthan=[N6index_Greaterthan i];
        end
        Tmax(i)=Tmax(i)+0;
    end

end
for i=1:num
    Tmax(i)=Tmax(i)+DeD(i);
end
  
    attack1=N1index(ceil(rand()*size(N1index,2)));
    attack2=N2index(ceil(rand()*size(N2index,2)));
    attack3=N3index(ceil(rand()*size(N3index,2)));
    attack4=N4index(ceil(rand()*size(N4index,2)));
    attack5=N5index(ceil(rand()*size(N5index,2)));
    attack6=N6index(ceil(rand()*size(N6index,2)));
    
    rn1=zeros(1,10);
    rn2=zeros(1,10);
    rn3=zeros(1,10);
    rn4=zeros(1,10);
    rn5=zeros(1,10);
    rn6=zeros(1,10);
    
    Norm1=zeros(1,10);
    Norm2=zeros(1,10);
    Norm3=zeros(1,10);
    Norm4=zeros(1,10);
    Norm5=zeros(1,10);
    Norm6=zeros(1,10);
    
    matrix1=matrix;
    matrix2=matrix;
    matrix3=matrix;
    matrix4=matrix;
    matrix5=matrix;
    matrix6=matrix;
    
    f1=DeD;
    f2=DeD;
    f3=DeD;
    f4=DeD;
    f5=DeD;
    f6=DeD;
    
    center1=[];center1=find(matrix1(attack1,:)==1);
    center2=[];center2=find(matrix2(attack2,:)==1);
    center3=[];center3=find(matrix3(attack3,:)==1); 
    center4=[];center4=find(matrix4(attack4,:)==1);
    center5=[];center5=find(matrix5(attack5,:)==1);
    center6=[];center6=find(matrix6(attack6,:)==1);
    
    centernum1=0;centernum1=size(find(matrix1(attack1,:)==1),2);
    centernum2=0;centernum2=size(find(matrix2(attack2,:)==1),2);
    centernum3=0;centernum3=size(find(matrix3(attack3,:)==1),2); 
    centernum4=0;centernum4=size(find(matrix4(attack4,:)==1),2);
    centernum5=0;centernum5=size(find(matrix5(attack5,:)==1),2);   
    centernum6=0;centernum6=size(find(matrix6(attack6,:)==1),2);
    
    matrix1(attack1,:)=0;matrix1(:,attack1)=0;
    matrix2(attack2,:)=0;matrix2(:,attack2)=0;
    matrix3(attack3,:)=0;matrix3(:,attack3)=0;
    matrix4(attack4,:)=0;matrix4(:,attack4)=0;
    matrix5(attack5,:)=0;matrix5(:,attack5)=0;
    matrix6(attack6,:)=0;matrix6(:,attack6)=0;
    A1=zeros(1,num);
    A2=zeros(1,num);
    A3=zeros(1,num);
    A4=zeros(1,num);
    A5=zeros(1,num);
    A6=zeros(1,num);
    for i=1:num
        A1(i)=sum(matrix1(i,:));
        A2(i)=sum(matrix2(i,:));
        A3(i)=sum(matrix3(i,:));
        A4(i)=sum(matrix4(i,:));
        A5(i)=sum(matrix5(i,:));   
        A6(i)=sum(matrix6(i,:));
    end

    ft1=f1;
    ft2=f2;
    ft3=f3;
    ft4=f4;    
    ft5=f5;
    ft6=f6;
    
    for i=1:centernum1
        ft1(1,center1(1,i))=f1(1,center1(1,i))+f1(1,attack1)/sum(A1(center1));
    end
    for i=1:centernum2
        ft2(1,center2(1,i))=f2(1,center2(1,i))+f2(1,attack2)/sum(A2(center2));
    end
    for i=1:centernum3
        ft3(1,center3(1,i))=f3(1,center3(1,i))+f3(1,attack3)/sum(A3(center3));
    end
    for i=1:centernum4
        ft4(1,center4(1,i))=f4(1,center4(1,i))+f4(1,attack4)/sum(A4(center4));
    end
    for i=1:centernum5
        ft5(1,center5(1,i))=f5(1,center5(1,i))+f5(1,attack5)/sum(A5(center5));
    end
    for i=1:centernum6
        ft6(1,center6(1,i))=f6(1,center6(1,i))+f6(1,attack6)/sum(A6(center6));
    end
    
    f1=ft1;
    f2=ft2;
    f3=ft3;
    f4=ft4;
    f5=ft5;
    f6=ft6;
    
    f1(1,attack1)=0;
    f2(1,attack2)=0;
    f3(1,attack3)=0; 
    f4(1,attack4)=0;
    f5(1,attack5)=0;
    f6(1,attack6)=0;
    
    Nb1=length(find(f1==0));
    Nb2=length(find(f2==0));
    Nb3=length(find(f3==0));
    Nb4=length(find(f4==0));
    Nb5=length(find(f5==0));
    Nb6=length(find(f6==0));
    rn1(1,1)=Nb1;
    rn2(1,1)=Nb2;
    rn3(1,1)=Nb3;
    rn4(1,1)=Nb4;
    rn5(1,1)=Nb5;
    rn6(1,1)=Nb6;
    
    [c1 d1]=components(sparse(matrix1));
    [c2 d2]=components(sparse(matrix2));
    [c3 d3]=components(sparse(matrix3));
    [c4 d4]=components(sparse(matrix4));
    [c5 d5]=components(sparse(matrix5));
    [c6 d6]=components(sparse(matrix6));
    Norm1(1,1)=max(d1);
    Norm2(1,1)=max(d2);
    Norm3(1,1)=max(d3);
    Norm4(1,1)=max(d4);
    Norm5(1,1)=max(d5);
    Norm6(1,1)=max(d6);
    
    te=0;
    while ~isempty(find(f1>Tmax(1,:)))
        attack1=find(f1>Tmax(1,:));
        center1=[];
        xnum1=zeros(size(attack1));
        for i=1:size(attack1,2) 
            center1=find(matrix1(attack1(1,i),:)==1);
            centernum1=size(center1,2);             
            matrix1(attack1(1,i),:)=0;
            matrix1(:,attack1(1,i))=0;            
            for k=1:num
                A1(k)=sum(matrix1(k,:));
            end
            ft1=f1;
            for j=1:centernum1
                ft1(1,center1(1,j))=f1(1,center1(1,j))+f1(1,attack1(1,i))/sum(A1(center1));
            end
            f1=ft1;
            f1(1,attack1(1,i))=0;       
            xnum1(i)=centernum1;
        end     
        Nb1=length(find(f1==0));
        te=te+1;
        if te==10
            break;
        end
        rn1(1,1+te)=Nb1;
        [c1 d1]=components(sparse(matrix1));
        Norm1(1,1+te)=max(d1);
    end
    
    te=0;
    while ~isempty(find(f2>Tmax(1,:)))
        attack2=find(f2>Tmax(1,:));
        center2=[];
        xnum2=zeros(size(attack2));
        for i=1:size(attack2,2) 
            center2=find(matrix2(attack2(1,i),:)==1);
            centernum2=size(center2,2);             
            matrix2(attack2(1,i),:)=0;
            matrix2(:,attack2(1,i))=0;            
            for k=1:num
                A2(k)=sum(matrix2(k,:));
            end
            ft2=f2;
            for j=1:centernum2
                ft2(1,center2(1,j))=f2(1,center2(1,j))+f2(1,attack2(1,i))/sum(A2(center2));
            end
            f2=ft2;
            f2(1,attack2(1,i))=0;       
            xnum2(i)=centernum2;
        end     
        Nb2=length(find(f2==0));
        te=te+1;
        if te==10
            break;
        end
        rn2(1,1+te)=Nb2;
        [c2 d2]=components(sparse(matrix2));
        Norm2(1,1+te)=max(d2);
    end
    
    te=0;
    while ~isempty(find(f3>Tmax(1,:)))
        attack3=find(f3>Tmax(1,:));
        center3=[];
        xnum3=zeros(size(attack3));
        for i=1:size(attack3,2) 
            center3=find(matrix3(attack3(1,i),:)==1);
            centernum3=size(center3,2);             
            matrix3(attack3(1,i),:)=0;
            matrix3(:,attack3(1,i))=0;            
            for k=1:num
                A3(k)=sum(matrix3(k,:));
            end
            ft3=f3;
            for j=1:centernum3
                ft3(1,center3(1,j))=f3(1,center3(1,j))+f3(1,attack3(1,i))/sum(A3(center3));
            end
            f3=ft3;
            f3(1,attack3(1,i))=0;       
            xnum3(i)=centernum3;
        end     
        Nb3=length(find(f3==0));
        te=te+1;
        if te==10
            break;
        end
        rn3(1,1+te)=Nb3;
        [c3 d3]=components(sparse(matrix3));
        Norm3(1,1+te)=max(d3);
    end
    
    te=0;
    while ~isempty(find(f4>Tmax(1,:)))
        attack4=find(f4>Tmax(1,:));
        center4=[];
        xnum4=zeros(size(attack4));
        for i=1:size(attack4,2) 
            center4=find(matrix4(attack4(1,i),:)==1);
            centernum4=size(center4,2);             
            matrix4(attack4(1,i),:)=0;
            matrix4(:,attack4(1,i))=0;            
            for k=1:num
                A4(k)=sum(matrix4(k,:));
            end
            ft4=f4;
            for j=1:centernum4
                ft4(1,center4(1,j))=f4(1,center4(1,j))+f4(1,attack4(1,i))/sum(A4(center4));
            end
            f4=ft4;
            f4(1,attack4(1,i))=0;       
            xnum4(i)=centernum4;
        end     
        Nb4=length(find(f4==0));
        te=te+1;
        if te==10
            break;
        end
        rn4(1,1+te)=Nb4;
        [c4 d4]=components(sparse(matrix4));
        Norm4(1,1+te)=max(d4);
    end
    
    te=0;
    while ~isempty(find(f5>Tmax(1,:)))
        attack5=find(f5>Tmax(1,:));
        center5=[];
        xnum5=zeros(size(attack5));
        for i=1:size(attack5,2) 
            center5=find(matrix5(attack5(1,i),:)==1);
            centernum5=size(center5,2);             
            matrix5(attack5(1,i),:)=0;
            matrix5(:,attack5(1,i))=0;            
            for k=1:num
                A5(k)=sum(matrix5(k,:));
            end
            ft5=f5;
            for j=1:centernum5
                ft5(1,center5(1,j))=f5(1,center5(1,j))+f5(1,attack5(1,i))/sum(A5(center5));
            end
            f5=ft5;
            f5(1,attack5(1,i))=0;       
            xnum5(i)=centernum5;
        end     
        Nb5=length(find(f5==0));
        te=te+1;
        if te==10
            break;
        end
        rn5(1,1+te)=Nb5;
        [c5 d5]=components(sparse(matrix5));
        Norm5(1,1+te)=max(d5);
    end
    
    te=0;
    while ~isempty(find(f6>Tmax(1,:)))
        attack6=find(f6>Tmax(1,:));
        center6=[];
        xnum6=zeros(size(attack6));
        for i=1:size(attack6,2) 
            center6=find(matrix6(attack6(1,i),:)==1);
            centernum6=size(center6,2);             
            matrix6(attack6(1,i),:)=0;
            matrix6(:,attack6(1,i))=0;            
            for k=1:num
                A6(k)=sum(matrix6(k,:));
            end
            ft6=f6;
            for j=1:centernum6
                ft6(1,center6(1,j))=f6(1,center6(1,j))+f6(1,attack6(1,i))/sum(A6(center6));
            end
            f6=ft6;
            f6(1,attack6(1,i))=0;       
            xnum6(i)=centernum6;
        end     
        Nb6=length(find(f6==0));
        te=te+1;
        if te==10
            break;
        end
        rn6(1,1+te)=Nb6;
        [c6 d6]=components(sparse(matrix6));
        Norm6(1,1+te)=max(d6);
    end
    
    LN1(y,:)=Norm1;
    LN2(y,:)=Norm2;
    LN3(y,:)=Norm3;
    LN4(y,:)=Norm4;
    LN5(y,:)=Norm5;
    LN6(y,:)=Norm6;
    
    NL1(y,:)=[0 rn1];
    NL2(y,:)=[0 rn2];
    NL3(y,:)=[0 rn3];
    NL4(y,:)=[0 rn4];
    NL5(y,:)=[0 rn5];
    NL6(y,:)=[0 rn6];
end

    NL1=NL1(:,1:10);
    NL2=NL2(:,1:10);
    NL3=NL3(:,1:10);
    NL4=NL4(:,1:10);
    NL5=NL5(:,1:10);
    NL6=NL6(:,1:10);
for y=1:100
   for i=1:10
       if NL1(y,i)~=0
           if NL1(y,i)~=36 && NL1(y,i)~=99
               NL1(y,i)=36-NL1(y,i);
           else
               NL1(y,i:end)=99;
           end
       end
       if NL2(y,i)~=0
           if NL2(y,i)~=36 && NL2(y,i)~=99
               NL2(y,i)=36-NL2(y,i);
           else
               NL2(y,i:end)=99;
           end
       end
       if NL3(y,i)~=0
           if NL3(y,i)~=36 && NL3(y,i)~=99
               NL3(y,i)=36-NL3(y,i);
           else
               NL3(y,i:end)=99;
           end
       end
       if NL4(y,i)~=0 
           if NL4(y,i)~=36 && NL4(y,i)~=99
               NL4(y,i)=36-NL4(y,i);
           else
               NL4(y,i:end)=99;
           end
       end
       if NL5(y,i)~=0
           if NL5(y,i)~=36 && NL5(y,i)~=99
               NL5(y,i)=36-NL5(y,i);
           else
               NL5(y,i:end)=99;
           end
       end
       if NL6(y,i)~=0
           if NL6(y,i)~=36 && NL6(y,i)~=99
               NL6(y,i)=36-NL6(y,i);
           else
               NL6(y,i:end)=99;
           end
       end
   end
    NL1(y,1)=36;
    NL2(y,1)=36;
    NL3(y,1)=36;
    NL4(y,1)=36;
    NL5(y,1)=36;
    NL6(y,1)=36;
end
for y=1:100
        if length(find(NL1(y,:)==0))>0
            k=length(find(NL1(y,:)==0));
            NL1(y,(10-k):end)=NL1(y,10-k);
        end
        if length(find(NL2(y,:)==0))>0
            k=length(find(NL2(y,:)==0));
            NL2(y,(10-k):end)=NL2(y,10-k);
        end
        if length(find(NL3(y,:)==0))>0
            k=length(find(NL3(y,:)==0));
            NL3(y,(10-k):end)=NL3(y,10-k);
        end
        if length(find(NL4(y,:)==0))>0
            k=length(find(NL4(y,:)==0));
            NL4(y,(10-k):end)=NL4(y,10-k);
        end
        if length(find(NL5(y,:)==0))>0
            k=length(find(NL5(y,:)==0));
            NL5(y,(10-k):end)=NL5(y,10-k);
        end
        if length(find(NL6(y,:)==0))>0
            k=length(find(NL6(y,:)==0));
            NL6(y,(10-k):end)=NL6(y,10-k);
        end
        for i=1:10
            if NL1(y,i)==99
                NL1(y,i)=0;
            end
            if NL2(y,i)==99
                NL2(y,i)=0;
            end
            if NL3(y,i)==99
                NL3(y,i)=0;
            end
            if NL4(y,i)==99
                NL4(y,i)=0;
            end
            if NL5(y,i)==99
                NL5(y,i)=0;
            end
            if NL6(y,i)==99
                NL6(y,i)=0;
            end
        end
end

NL1_per=sum(NL1)/(num*100);
NL2_per=sum(NL2)/(num*100);
NL3_per=sum(NL3)/(num*100);
NL4_per=sum(NL4)/(num*100);
NL5_per=sum(NL5)/(num*100);
NL6_per=sum(NL6)/(num*100);

figure;
plot(1:10,NL1_per,'-k+',1:10,NL2_per,'-g*',1:10,NL3_per,'-b^',1:10,NL4_per,'-rh',1:10,NL5_per,'-cs',1:10,NL6_per,'-mo','linewidth',2,'MarkerSize',10);
axis([1 10 0 1]);
set(gca, 'FontName', 'Times New Roman','Fontsize', 10.5);
legend('k=1','k=2','k=3','k=4','k=5','k=6');
xlabel('Iteration');ylabel('Node Livability');