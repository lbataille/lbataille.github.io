function C = corr2neigh(A,B,se)
    if (ismember(numel(size(A)),[2;3]) && ismember(numel(size(A)),[2;3]))
        if ((size(A,1)==size(B,1)) && (size(A,2)==size(B,2)))
            if xor(numel(size(A))==3,numel(size(B))==3)
                C = nan(size(A));        
                if (numel(size(A)) ==3)

                    for i = 1:size(A,3)
                        C(i,:) = corr2neigh(A(:,:,i),B,se);
                    end
                elseif (numel(size(B)) ==3)
                                C = nan(size(A));
                    for i = 1:size(A,3)
                        C(i,:) = corr2neigh(A,B(:,:,i),se);
                    end
                end
            elseif ((numel(size(A))==3) && (numel(size(B))==3))
                C = nan(size(A));        
                for i = 1:size(A,3)
                    C(i,:) = corr2neigh(A(:,:,i),B(:,:,i),se);
                end     
            elseif ((numel(size(A))==2) && (numel(size(B))==2)) 
                h0 = se.Neighborhood;
                h1 = h0/nnz(h0);

                Am = filter2(h1,A);
                Bm = filter2(h1,B);
                Astd = (filter2(h1,(A-Am).^2)).^0.5;
                Bstd = (filter2(h1,(B-Bm).^2)).^0.5;

                ABcov = filter2(h1,(A-Am).*(B-Bm));
                C = ABcov./(Astd.*Bstd);
                C((Astd==0) & (Bstd==0)) = 1;
                C(xor((Astd==0),(Bstd==0))) = 0;
            end
        end
    end
    C(~isnan(C) & abs(C)>1) = nan;
end