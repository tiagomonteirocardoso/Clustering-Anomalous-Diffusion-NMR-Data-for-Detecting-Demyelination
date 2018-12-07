function outmatrix = burning_algo( matrix,resol )

connectmatrix = matrix;
cod=1000;
res=resol;

for p1=1:res
    for q1=1:res
        for r1=1:res
            if connectmatrix(p1,q1,r1)==0
                fila=[p1,q1,r1];
                cod=cod+1;
            
                while numel(fila) ~= 0
                    fila = unique(fila,'rows');
                    label(fila);
                    fila = neighbours(fila);
                end
            end
        end
    end
end

    function label(ctl)
        
        [numrows numcols] = size(ctl);
        
         for k=1:numrows
             connectmatrix(ctl(k,1),ctl(k,2),ctl(k,3))=cod;
         end
        
    end

    function outvec=neighbours(invec)
        
        [numrows2 numcols2] = size(invec);
        indices=[];
        count=1;
        
        for l=1:numrows2
            a=invec(l,1);
            b=invec(l,2);
            c=invec(l,3);
            
            if a==1
                nr=[a,a+1];
            elseif a==res
                nr=[a-1,a];
            else
                nr=[a-1,a,a+1];
            end
            
            if b==1
                nc=[b,b+1];
            elseif b==res
                nc=[b-1,b];
            else
                nc=[b-1,b,b+1];
            end
            
            if c==1
                nh=[c,c+1];
            elseif c==res
                nh=[c-1,c];
            else
                nh=[c-1,c,c+1];
            end
            
            num_rows=length(nr);
            num_cols=length(nc);
            num_floors=length(nh);
            
            for m=1:num_rows
                for n=1:num_cols
                    for o=1:num_floors
                        if ((nr(m)==a)+(nc(n)==b)+(nh(o)==c)==2)&(connectmatrix(nr(m),nc(n),nh(o))==0)
                            indices(count,:) = [nr(m),nc(n),nh(o)];
                            count=count+1;
                        end
                    end
                end
            end
        end
        
        outvec=indices;
        indices=[];
        
    end

outmatrix = connectmatrix;
end