function [thr,nClass, prob] = isBimodMatr(B,prob_cut,diff_cut)
        
        n_curr = size(B,1);

        Bvec = sort(B(~eye(n_curr)));
        mvec = min(Bvec);
        Blog = log(B);
        Blog = sort(Blog(~eye(n_curr)));
        Bdiff = diff(Blog);
        max_diff = max(Bdiff);
        mind = find(Bdiff==max_diff,1,'first');
        thr = (Bvec(mind)+Bvec(mind+1))/2;
        
        ind1 = Bvec < thr;
        ind2 = Bvec >= thr;
        if sum(ind1) >= sum(ind2)
            try
                phat = betafit(Bvec(ind1));
            catch
                phat = logbetafit(Bvec(ind1));
            end
            prob = 1 - betacdf(Bvec(find(ind2,1,"first")),phat(1),phat(2));
        else
            try
                phat = betafit(Bvec(ind2));
            catch
                phat = logbetafit(Bvec(ind2));
            end
            prob = betacdf(Bvec(find(ind1,1,"last")),phat(1),phat(2));
        end
        if mvec > 1e-7
            nClass = 2;
            return;
        end
        if (prob > prob_cut)||(sum(ind1)<0.2*length(ind1)) %||(max_diff < diff_cut)
            nClass = 1;
        else
            nClass = 2;
        end
