tic
updates=100;
itt = 100000;
update_frac=itt/updates;
N_eff=round(itt/update_frac);
parfor_progress_imp(updates);
parfor i=1:itt
    pause(10*rand/itt); % Replace with real code
    if rand<(1/update_frac)
        parfor_progress_imp;
    end
end
parfor_progress_imp(0);
toc
