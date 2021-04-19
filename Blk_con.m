classdef Blk_con <handle
properties
    defName
    alias
    database
    hash

    srcTable
    srcKey

    ptchsOpts
    blkOpts

    %%
    cmpLookup
    cmpKey

    cndLookup
    cndKey

    lvlLookup
    lvlKey
    %

    optsTable
    optsKey

    blkTable
    blkKey

    selInd

    %
    dims
    nDim
    nCmpPerDim
    nLvlPerDim
end
properties(Hidden)
    stdRC
    cmpRC
    replaceBind
end
methods(Access = ?Blk)
    function obj=Blk_con(defName,replaceBind,bSave)
        if ~startsWith(defName,'D_exp_')
            defName=['D_blk_' defName];
        end
        if exist('replaceBind','var') && ~isempty(replaceBind)
            obj.replaceBind=replaceBind;
        end
        if ~exist('bSave','var') || isempty(bSave)
            bSave=1;
        end
        fname=which(defName);
        if isempty(fname)
            error('Cant find file from defName');
        end
        obj.defName=defName;
        obj.alias=sed('s',defName,'^D_blk_','');

        [blkOpts,ptchsOpts]=obj.load_def_();
        obj.load_src_table_();

        obj.parse_blkOpts_(blkOpts);
        obj.parse_ptchsOpts_(ptchsOpts);

        obj.get_opts_tables_();
        obj.get_lvlInd_tables_();
        obj.get_blk_table_();
        if bSave && ~isempty(obj.replaceBind)
            obj.move();
        end
        obj.save();
    end
    function obj=move(obj)
        dire=Blk.get_dir(obj.alias);
        table=obj.optsTable;
        key=obj.optsKey;

        n=[];
        while true

            foptso=[dire 'opts' '.mat'];
            fopts=[dire 'opts_old' num2str(n) '.mat '];

            flko=[dire 'lookup' '.mat'];
            flk=[dire 'lookup_old' num2str(n) '.mat '];

            fblko=[dire 'blk' '.mat'];
            fblk=[dire 'blk_old' num2str(n) '.mat '];
            if ~exist(fopts,'file') && ~exist(flk,'file') && ~exist(fblk,'file')
                movefile(foptso,fopts);
                movefile(flko,flk);
                movefile(fblko,fblk);
                break
            elseif isempty(n)
                n=1;
            else
                n=n+1;
            end
        end

    end
    function obj=save(obj)
        dire=Blk.get_dir(obj.alias);
        if ~exist(dire,'dir')
            mkdir(dire);
        end

        table=obj.optsTable;
        key=obj.optsKey;
        save([dire 'opts'],'table','key');

        table=obj.blkTable;
        key=obj.blkKey;
        save([dire 'blk'],'table','key');

        lookup=Blk.make_lookup_struct(obj);
        save([dire 'lookup'],'lookup');


    end
end
methods(Access=protected)
%% UTIL
    function obj=load_src_table_(obj)
        [obj.srcTable,obj.srcKey]=imapPch.load_src_table(obj.database,obj.hash);
    end
    function col=get_src_column(obj,key)
        ind=ismember(obj.srcKey,key);
        col=obj.srcTable(:,ind);
    end
    function col=get_block_column(obj,key)
        if ~iscell(key)
            key={key};
        end
        gd=ismember(key,obj.blkKey);
        if all(gd)
            ind=ismember(obj.blkKey,key);
            col=obj.blkTable(:,ind);
            return
        end

        col=zeros(size(obj.blkTable,1), numel(key));
        ind=ismember(obj.blkKey,key(gd));
        col(:,gd)=obj.blkTable(:,ind);
        col(:,~gd)=obj.get_block_column_dim_lvl(key(~gd));
    end
    function col=get_block_column_dim_lvl(obj,name)
        col=zeros(size(obj.blkTable,1), numel(name));
        for i = 1:length(name)
            d=find(ismember(obj.dims,name{i}));
            lvls=obj.get_block_column('lvlInd');
            a=obj.lvlInd_to_lvlRC(lvls);
            col(:,i)=a(:,d);
        end
    end
    function col=get_block_column_dim_cmp(obj,name)
        if ~iscell(name)
            name={name};
        end
        col=zeros(size(obj.blkTable,1), numel(name));
        for i = 1:length(name)
            d=find(ismember(obj.dims,name{i}));
            lvls=obj.get_block_column('lvlInd');
            a=obj.lvlInd_to_lvlRC(lvls);
            col(:,i)=a(:,d);
        end
    end
    function cnd=lvlInd_cmpInd_to_cndInd(obj,lvl,cmp)
        nStd=size(obj.lvlLookup,1);
        nCmp=size(obj.cmpLookup,1);

        cnd=sub2ind([nCmp nStd],cmp,lvl);

    end
    function lvl=cndInd_to_lvlInd(obj,cnd)
        lvl=[obj.cndLookup(cind,2)];
    end
    function cmps=cndInd_to_cmpInd(obj,cnd)
        cmps=[obj.cndLookup(cind,3)];
    end

    %% IND2RC
    function RC=lvlInd_to_lvlRC(obj,ind)
        RC=[obj.lvlLookup(ind,2:end)];
    end
    function RC=cmpInd_to_cmpRC(obj,ind)
        RC=[obj.cmpLookup(ind,2:end)];
    end
%% PARSE
    function obj=parse_blkOpts_(obj,blkOpts)
        obj.blkOpts=Blk_con.parse_blkOpts(blkOpts);
    end
    function obj=parse_ptchsOpts_(obj,ptchsOpts)
        obj.ptchsOpts=Blk_con.parse_ptchsOpts(ptchsOpts);
        obj.dims=obj.ptchsOpts.dims;
    end
    function [blkOpts,ptchsOpts]=load_def_(obj,def)
        run(obj.defName);

        if ~exist('dtb','var')
            error('Database name required in def-file');
        end
        obj.database=dtb;

        if exist('pchAlias','var') && ~isempty(pchAlias)
            hashes=imapCommon.alias2hashes(pchAlias,obj.database);
            obj.hash=hashes.pch;
        elseif exist('hash','var') && ~isempty(hash)
            obj.hash=hash;
        else
            error('Patches hash or alias required in def-file');
        end

        blkOpts=struct();
        flds=Blk_con.get_blk_flds();
        for i = 1:length(flds)
            fld=flds{i};
            if exist(fld,'var')
                blkOpts.(fld)=eval(fld);
            end
        end

        ptchsOpts=struct();
        flds=Blk_con.get_ptchs_flds();
        for i = 1:length(flds)
            fld=flds{i};
            if exist(fld,'var')
                ptchsOpts.(fld)=eval(fld);
            end
        end
    end
%% Tables
    function obj=get_opts_tables_(obj)
        [obj.optsTable,obj.optsKey, obj.nDim,obj.nLvlPerDim,obj.nCmpPerDim, obj.stdRC,obj.cmpRC]=Blk_con.get_opts_tables(obj.ptchsOpts);


    end
    function obj=get_lvlInd_tables_(obj)
        [~,~,stdInd]=unique(obj.stdRC,'rows');
        [~,~,cmpInd]=unique(obj.cmpRC,'rows');

        obj.lvlLookup=unique([stdInd obj.stdRC],'rows');
        obj.cmpLookup=unique([cmpInd obj.cmpRC],'rows');

        obj.lvlKey=['stdInd',obj.dims];
        obj.cmpKey=['cmpInd',obj.dims];

        obj.cndLookup=[transpose(1:length(stdInd)),stdInd,cmpInd];
        obj.cndKey={'cndInd','lvlInd','cmpInd'};


    end
    function obj=get_blk_table_(obj)
        b=obj.blkOpts;
        [obj.blkTable,obj.blkKey]=Blk_con.get_blk_table(obj.nDim,obj.nLvlPerDim,obj.nCmpPerDim,b.modes,b.nBlkPerDim,b.nTrlPerLvl,b.nIntrvlPerTrl,b.sd);

        %get lvl and cmp ind
        I=find(ismember(obj.blkKey,'lvlInd'));
        lvlInd=obj.blkTable(:,I);
        cmpInd=obj.get_cmpInd_blk_();

        % remove stdind
        obj.blkTable(:,I)=[];
        obj.blkKey(I)=[];

        cndInd=obj.lvlInd_cmpInd_to_cndInd(lvlInd,cmpInd);

        %cndInd=zeros(size(lvlInd));
        %for i = 1:size(lvlInd,1)
        %    cndInd(i)=find(ismember(cndTable,old(i,:),'rows'));
        %end

        obj.blkTable=[cndInd lvlInd cmpInd obj.blkTable];
        obj.blkKey=['cndInd' 'lvlInd' 'cmpInd' obj.blkKey];


        % APPEND P
        obj.get_selInd_table_blk_([]);
        if ~isempty(obj.replaceBind)
            obj.get_selInd_table_blk_(obj.replaceBind);
        end
        obj.blkTable=[obj.blkTable obj.selInd];
        obj.blkKey=[obj.blkKey 'P'];

    end
    function cmpInd=get_cmpInd_blk_(obj)


        b=obj.blkOpts;
        nModes=numel(b.modes);
        cmpInd=Blk_con.get_cmpInd_blk(nModes,obj.nLvlPerDim,obj.nCmpPerDim,b.nBlkPerDim,b.nTrlPerLvl,b.nIntrvlPerTrl);
    end
    function obj=get_selInd_table_blk_(obj,replaceBind)
        %Blk key={'mode','lvlInd','blk','trl','intrvl','cmpInd','cmpNum'};

        b=obj.blkOpts;
        repeats=b.repeats;
        mirror=b.mirror;
        pmirror=b.pmirror;
        binVals=obj.get_src_column('B');
        bins=[obj.ptchsOpts.bins{:}];
        binVals=vertcat(binVals{:});
        if exist('replaceBind','var') && ~isempty(replaceBind)
            idx=obj.get_src_column('P');
            idx=vertcat(idx{:});
            size(replaceBind)
            size(idx)
            binGd=~replaceBind & ismember(idx,obj.selInd); % sample from
            replaceBind=ismember(obj.selInd, idx(replaceBind));
        else
            binGd=true(size(binVals));
        end

        nStd=prod(obj.nLvlPerDim);
        nModes=numel(b.modes);
        nBins=numel(bins);
        [nIntrvlPerBin,nIntrvlAll]=Blk_con.get_nIntrvlPerBin(nBins,nModes,obj.nLvlPerDim,obj.nCmpPerDim,b.nBlkPerDim,b.nTrlPerLvl,b.nIntrvlPerTrl);

        % XXX
        binCol=obj.get_block_column('bins');
        %nIntrvlPerBin 27000
        mr=[mirror pmirror repeats];
        for i = 1:length(mr)
            m=mr{i};
            switch m
            case 'lvlInd'
                nIntrvlPerBin=nIntrvlPerBin/nStd;
            case obj.dims
                ind=ismember(obj.dims,m);
                nIntrvlPerBin=nIntrvlPerBin/obj.nLvlPerDim(ind);
            case 'mode'
                nIntrvlPerBin=nIntrvlPerBin/nModes;
            otherwise obj.dims
                error(['Unhandled repeat or mirror case: ' m ]);
            end
        end
        disp(['nIntrvlPerBin ' num2str(nIntrvlPerBin)]);
        % nIntrvlPerBin1 1800

        if ~exist('replaceBind','var') || isempty(replaceBind)
            replaceBind=true(size(binCol));
            obj.selInd=zeros(nIntrvlAll,1);
            bReplaceFlag=0;
        else
            bReplaceFlag=1;
        end


        for i = 1:length(bins)
            binBind=binCol==i & replaceBind; % Available to set
            obj=bin_table_fun(obj,bins(i),binVals,binBind,repeats,mirror,pmirror,nIntrvlPerBin,binGd,bReplaceFlag);
        end

        function obj=bin_table_fun(obj,bin,binVals,binBind,repeats,mirror,pmirror,nIntrvlPerBin,binGd,bReplaceFlag)

            % TODO
            % bins
            % and/or?


            binInd=find(ismember(binVals,bin) & binGd); %Available to SAMPLE FROM
            nStm=numel(binInd); % 7008

            bReplace=false;
            if nStm < nIntrvlPerBin & ~bReplaceFlag
                bReplace=true;
                disp([ 'Setting replace to true. N aviablable stim in bin ' num2str(bin) ' = ' num2str(nStm) '. Req = ' num2str(nIntrvlPerBin) ]);
            end

            if isempty(repeats)
                nUnqRep=1;
                uIndsRep=true(size(binBind));
            else
                col=obj.get_block_column(repeats);
                [~,~,uIndsRep]=unique(col,'rows');
                nUnqRep=(max(uIndsRep));
            end

            if isempty(mirror)
                nUnqMir=1;
                uIndsMir=true(size(binBind));
            else
                col=obj.get_block_column(mirror);
                [u,~,uIndsMir]=unique(col,'rows');
                nUnqMir=numel(u);
            end

            if isempty(pmirror)
                nUnqPmir=1;
                uIndsPmir=true(size(binBind));
            else
                col=obj.get_block_column(pmirror);
                [u,~,uIndsPmir]=unique(col,'rows');
                nUnqPmir=numel(u);
            end

            %% mirror
            indsMir=false(size(binBind));
            %75000/3 modes /5 bins / 5stds
            for i = 1:nUnqRep
                ind=(uIndsRep==i) & (uIndsMir==1) & (uIndsPmir==1) & binBind;
                if sum(ind)==0
                    continue
                end
                obj.selInd(ind)=datasample(binInd,sum(ind),'Replace',bReplace);
            end

            indsMir=uIndsMir & binBind;
            for i = 2:nUnqMir
                ind=(uIndsMir==i) & binBind;
                obj.selInd(ind)=obj.selInd(indsMir);
            end

            for j = 1:nUnqRep
                indsPmir=(uIndsRep==j) & (uIndsPmir==1) & binBind;
                K=sum(indsPmir);
                data=obj.selInd(indsPmir);
                for i = 2:nUnqPmir
                    ind=(uIndsRep==j) & (uIndsPmir==i) & binBind;
                    obj.selInd(ind)=data(randperm(K));
                end
            end

        end
    end
end
methods(Static,Access=protected)
    function B=shuffle_within_rows_(A)
        [M,N] = size(A);
        % Preserve the row indices
        rowIndex = repmat(transpose(1:M),[1 N]);
        % Get randomized column indices by sorting a second random array
        [~,randomizedColIndex] = sort(rand(M,N),2);
        % Need to use linear indexing to create B
        newLinearIndex = sub2ind([M,N],rowIndex,randomizedColIndex);
        B = A(newLinearIndex);
    end
    function cmpInd=get_cmpInd_blk(nModes,nLvlPerDim,nCmpPerDim,nBlkPerDim,nTrlPerLvl,nIntrvlPerTrl)
        nStd=prod(nLvlPerDim);
        nCmp=prod(nCmpPerDim);
        nBlk=nModes.*nStd.*nBlkPerDim; % 5 * 5 * 3
        nTrl=nTrlPerLvl.*nCmp;
        nTrlPerBlk=nTrl./nBlkPerDim; % 180
        cmpSubs=Blk_con.get_cmpSubs(nCmpPerDim,nTrlPerBlk,nBlk,nIntrvlPerTrl);
        [~,~,cmpInd]=unique(cmpSubs,'rows');
    end
    function [table,key]=get_blk_table(nDim,nLvlPerDim,nCmpPerDim,modes,nBlkPerDim,nTrlPerLvl,nIntrvlPerTrl,sd)

        nStd=prod(nLvlPerDim);
        nCmp=prod(nCmpPerDim);
        nModes=numel(modes);
        nBlk=nModes.*nStd.*nBlkPerDim; %50
        nTrl=nTrlPerLvl.*nCmp;
        nTrlPerBlk=nTrl./nBlkPerDim; % 180

        % 3 25 5 100 2 XXX
        %nTrlPerBlk
        %nBlkPerDim
        %nIntrvlPerTrl
        %nModes*nStd*nBlkPerDim*nTrlPerBlk*nIntrvlPerTrl
        % 1 point = 100 trials
        % 1 curve = 500 trials = 1000 uique stim/bin
        % across each bin = 1000*5 = 5000 unique stim total
        % across each dsp = 5000*5 = 25000 stim total
        % 75000 for all modes
        % 15000 each bin
        table=distribute(1:nModes,1:nStd,1:nBlkPerDim,1:nTrlPerBlk,1:nIntrvlPerTrl); % 18000

        rng(sd);
        cmpNum=Blk_con.get_cmp_num(nTrlPerBlk,nBlk,nIntrvlPerTrl);

        table=[table cmpNum];
        key={'mode','lvlInd','blk','trl','intrvl','cmpNum'};
    end
    function c=get_cmpSubs(nCmpPerDim,nTrlPerBlk,nBlk,nIntrvlPerTrl)
        % which comparison to use. 1-5 eg ncol = n dims
        nIntrvlAll=nTrlPerBlk*nBlk*nIntrvlPerTrl;
        c= zeros(nIntrvlAll,length(nCmpPerDim));
        for i = 1:length(nCmpPerDim)
            c(:,i)=cmp_fun(nCmpPerDim(i),nTrlPerBlk,nBlk,nIntrvlPerTrl);
        end
        function c=cmp_fun(nCmpPerDim,nTrlPerBlk,nBlk,nIntrvlPerTrl)
            c=repmat(repelem(1:nCmpPerDim,1,nTrlPerBlk/nCmpPerDim),nBlk,1);
            c=transpose(Blk_con.shuffle_within_rows_(c));
            c=c(:);
            counts=hist(c(1:nTrlPerBlk),unique(c(1:nTrlPerBlk)));
            if ~isuniform(counts)
                error('something bad happend');
            end
            c=repelem(c,nIntrvlPerTrl,1);
        end
    end
    function c=get_cmp_num(nTrlPerBlk,nBlk,nIntrvlPerTrl)
        % std or comparison 1 or comparison 2 ...

        c=repmat((1:nIntrvlPerTrl),nTrlPerBlk*nBlk,1);
        c=transpose(Blk_con.shuffle_within_rows_(c));
        c=c(:)-1;
    end
end
methods(Static)
    function opts=parse_ptchsOpts(opts)
        % XXX TODO fill in missing params if not exist
        %opts=parse([],Opts,P);

        P=Blk_con.get_ptchs_parseOpts();
        flds=P(:,1);

        opts=parse_inds_fun(opts,P);

        function opts=parse_inds_fun(opts,P)
            flds=P(:,1);
            Sz=zeros(size(P,1),1);
            for i = 1:length(flds)
                fld=flds{i};
                if ~isfield(opts,fld)
                    opts.(fld)=P{i,2};
                end
            end
            for i = 1:size(P,1)
                fld=flds{i};
                val=opts.(fld);

                if iscell(val)
                    for j = 1:numel(val)
                        opts.(fld){j}=parse_ind_fun(val{j},fld);
                    end
                else
                    opts.(fld)=parse_ind_fun(val,fld);
                end
            end
        end
        function val=parse_ind_fun(val,fld)
            sz=size(val);
            if isempty(val);
                return
            end

            if endsWith(fld,'PosXYZm')
                if sz(1) == 3 && sz(2) ~= 3
                    val=transpose(val);
                elseif sz(1) ~=3 && sz(2) ~=3
                    error([ fld ' second dimension must be size 3']);
                end
            elseif endsWith(fld,'XYdeg') || strcmp(fld,'wdwSymInd')
                if sz(1) == 2 && sz(2) ~= 2
                    val=transpose(val);
                elseif sz(1) ~=2 && sz(2) ~=2
                    error([ fld ' second dimension must be size 2']);
                end
            elseif sz(1) >= 1 && sz(2) == 1
                   ;
            elseif ~ischar(val) && sz(1) == 1 && sz(2) > 1
                val=transpose(val);
            elseif ~ischar(val)
                error([ fld ' second dimension must be size 1']);
            end

            if iscell(val) && ~all(cellfun(@(x) ismember(x,{'disp','win'}),val));
                error([ fld ' can only be ''disp'' or ''win''.']);
            elseif ischar(val) && endsWith(fld,'DispOrWin') && ~ismember(val,{'disp','win'})
                error([ fld ' can only be ''disp'' or ''win''.']);
            elseif  endsWith(fld,'DispOrWin') && ~ischar(val) && ~iscell(val)
                error([ fld ' can only be ''disp'' or ''win''.']);
            end
        end
    end

    function [table,key,nRow,nCol,nAisle, stdRC,cmpRC]=get_opts_tables(opts)
        P=Blk_con.get_ptchs_parseOpts();
        flds=P(:,1);

        [nRow,nCol,nAisle]=get_dimensions_fun(opts,flds);

        [stdRC,cmpRC]=get_RC(nCol,nAisle);

        [table,key]=expand_fun(opts,flds,stdRC,cmpRC);


        function [table,key]=expand_fun(opts,flds,stdRC,cmpRC)
            dimflds=flds(ismember(flds,opts.dims));
            flds(ismember(flds,[{'linked','dims'} opts.dims]))=[];
            N=size(stdRC,1);
            M=length(flds);
            ndim=size(stdRC,2);

            table=cell(N, M);
            for j = 1:M
                fld=flds{j};
                if ismember(fld,horzcat(opts.linked{:}))
                    table(:,j)=expand_cell_fun(opts,fld);
                elseif iscell(opts.(fld))
                    error('Cells are reserved for linked params');
                else
                    table(:,j)=expand_single_fun(opts,fld,N);
                end
            end

            dimTable=cell(N,ndim);
            for i = 1:length(dimflds)
                fld=dimflds{i};

                ind=find(ismember(opts.dims,fld));
                [stds cmps] = expand_dim_fun(opts,fld,ind,stdRC,cmpRC,N);

                if size(stds,1) == 1 & size(stds,2) > 1
                    stds=transpose(stds);
                end
                if size(cmps,1) == 1 & size(cmps,2) > 1
                    cmps=transpose(cmps);
                end

                dim=mat2cell([stds cmps],ones(size(stds,1),1));
                dimTable(:,ind)=dim;
            end

            table=[dimTable table];
            f=transpose(flds);
            f(ismember(f,opts.dims))=[];
            key=[opts.dims f];
        end
        function [std,cmp]=expand_dim_fun(opts,fld,ind,stdRC,cmpRC,N);
            c=stdRC(:,ind);
            a=cmpRC(:,ind);
            if max(a) > 1
                a=a+1;
            end
            F=size(opts.(fld));
            o=ones(size(c));

            ind=sub2ind(F, o,c,o);
            std=opts.(fld)(ind);

            ind=sub2ind(F, o,c,a);
            cmp=opts.(fld)(ind);

        end
        function val=expand_single_fun(opts,fld,N)
            val=repmat({opts.(fld)},N,1);
        end

        function [stdRC,cmpRC]=get_RC(nCol,nAisle)
        %std{:} cmp{:}
            C=cellfun(@(x) 1:x, num2cell(nCol),UO,false);
            % every standard combination

            A=cellfun(@(x) 1:x, num2cell(nAisle),UO,false);
            % every cmp combination
            if numel(C) == 1
                C=C{1};
            else
                C=distribute(C{:});
            end
            if numel(A) == 1
                A=A{1};
            else
                A=distribute(A{:});
            end

            RC=distribute(C,A);
            % every combination of standard combinations and cmp combinations

            h=size(RC,2)/2;
            stdRC=RC(:,1:h);
            cmpRC=RC(:,h+1:end);
        end
        function [nRow,nCol,nAisle]=get_dimensions_fun(opts,flds)
            nRow=numel(opts.dims);
            nCol=zeros(nRow,1);
            nAisle=zeros(nRow,1);
            for i = 1:length(opts.dims)
                nCol(i)=size(opts.(opts.dims{i}),2);
                nAisle(i)=size(opts.(opts.dims{i}),3);
            end
            ind=nAisle > 1;
            nAisle(ind)=nAisle(ind)-1;
        end
    end
    function flds=get_blk_flds();
        P=Blk_con.get_blk_parseOpts();
        flds=P(:,1);
    end
    function flds=get_ptchs_flds();
        P=Blk_con.get_ptchs_parseOpts();
        flds=P(:,1);
    end
    function P=get_blk_parseOpts()
        P={...
              'dims',[],[] ...
              ;'repeats',[],'' ...
              ;'mirror',[],'' ...
              ;'pmirror',[],'' ...
              ;'modes',[],'isallint' ...
              ;'nBlkPerDim',[],'isint' ...
              ;'nTrlPerLvl',[],'isint' ...
              ;'nIntrvlPerTrl',[],'isint' ...
              ;'sd',[],'isint', ...
          };
    end
    function P=get_ptchs_parseOpts()
        P={...
             'dims',{},'iscell_e' ...
            ;'linked',{},'iscell_e' ...
            ;'disparity',[],'isallnum_e'  ... % 1
            ;'speed',[],'isallnum_e'  ... % 1
            ;'bins',[],'isallint_e' ...
            ...
            ;'rmsFix',[],'isallnum_e'  ... % 1
            ;'dcFix',[],'isallnum_e'  ... % 1
            ;'dnkFix',[],'isallnum_e'  ... % 1
            ...
            ;'trgtDispOrWin','disp',''  ... %1, 'disp' 'win'
            ;'trgtPosXYZm',[0 0 0],'isallnum_e'  ...   % 3
            ...
            ;'focDispOrWin','disp',''  ... %1, 'disp' 'win'
            ;'focPosXYZm',[0 0 0],'isallnum_e'  ...   % 3
            ...
            ;'stmPosXYZm',[0 0 0],'isallnum_e'  ... %1, 'disp' 'win'
            ;'stmXYdeg',[],'isallnumeric_e'  ...   % 3
            ...
            ;'wdwType',[],'ischar'  ...   %1
            ;'wdwPszRCT',[],''  ...   % 2-3 OR char
            ;'wdwPszRCT',[],''  ...   % 2-3 OR char
            ;'wdwRmpDm',[],''  ...    % 1-3
            ;'wdwDskDm',[],''  ...    % 1-3
            ;'wdwSymInd',[],''  ...   % 1-3
             ...
            ;'duration',[],''  ...   % 1-3

          };
    end
    function out=parse_blkOpts(Opts)
        P=Blk_con.get_blk_parseOpts();
        out=parse([],Opts,P);
    end
    function [nIntrvlPerBin,nIntrvlAll]=get_nIntrvlPerBin(nBins,nModes,nLvlPerDim,nCmpPerDim,nBlkPerDim,nTrlPerLvl,nIntrvlPerTrl)
        nStd=prod(nLvlPerDim);
        nIntrvlAll=nTrlPerLvl*prod(nCmpPerDim)*nStd*nIntrvlPerTrl*nModes;
        nIntrvlPerBin=nIntrvlAll/nBins;
    end
end
end
