classdef Blk <handle
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

    lvlIndTable
    lvlIndKey

    optsTable
    optsKey

    blkTable
    blkKey

    selInd

    dims
    nDim
    nCmpPerLvl
    nLvlPerDim
end
properties(Hidden)
end
methods
    function obj=Blk(defName)
        if ~startsWith(defName,'D_exp_')
            defName=['D_blk_' defName];
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

        obj.get_blk_table_();
        obj.get_selInd_table_();
        size(obj.blkTable)
        size(obj.optsTable)
        size(obj.selInd)
        size(obj.lvlIndTable)
        obj.save();
    end
    function obj=save(obj)
        % XXX
    end
    function col=get_block_column(obj,key)
        ind=ismember(obj.blkKey,key);
        col=obj.blkTable(:,ind);
    end
    function col=get_src_column(obj,key)
        ind=ismember(obj.srcKey,key);
        col=obj.srcTable(:,ind);
    end
end
methods(Access=private)
    function obj=parse_blkOpts_(obj,blkOpts)
        obj.blkOpts=Blk.parse_blkOpts(blkOpts);
    end
    function obj=parse_ptchsOpts_(obj,ptchsOpts)
        obj.ptchsOpts=Blk.parse_ptchsOpts(ptchsOpts);
        obj.dims=obj.ptchsOpts.dims;
    end
    function [blkOpts,ptchsOpts]=load_def_(obj,def)
        run(obj.defName);

        if ~exist('dtb','var')
            error('Database name required in def-file');
        end
        obj.database=dtb;

        if ~exist('hash','var') || isempty(hash)
            error('Patches hash required in def-file');
        end
        obj.hash=hash;

        blkOpts=struct();
        flds=Blk.get_blk_flds();
        for i = 1:length(flds)
            fld=flds{i};
            if exist(fld,'var')
                blkOpts.(fld)=eval(fld);
            end
        end

        ptchsOpts=struct();
        flds=Blk.get_ptchs_flds();
        for i = 1:length(flds)
            fld=flds{i};
            if exist(fld,'var')
                ptchsOpts.(fld)=eval(fld);
            end
        end
    end
    function obj=load_src_table_(obj)
        [obj.srcTable,obj.srcKey]=imapPch.load_src_table(obj.database,obj.hash);
    end
%% Tables
    function obj=get_opts_tables_(obj)
        [obj.optsTable,obj.optsKey, obj.nDim,obj.nLvlPerDim,obj.nCmpPerLvl, stdRC,cmpRC,cndRC]=Blk.get_opts_tables(obj.ptchsOpts);

        obj.lvlIndTable=[stdRC,cmpRC,cndRC];
        obj.lvlIndKey=[ ...
                      'subInd','lvlInd' ...
                      repmat({'stdInd'},1,size(stdRC,2)) ...
                      repmat({'cmpInd'},1,size(cmpRC,2)) ...
                    ];...
    end
    function obj=get_blk_table_(obj)
        b=obj.blkOpts;
        [obj.blkTable,obj.blkKey]=Blk.get_blk_table(obj.nDim,obj.nLvlPerDim,obj.nCmpPerLvl,b.modes,b.nBlkPerLvl,b.nTrlPerLvl,b.nIntrvlPerTrl,b.sd);
        obj.get_selInd_table_();
    end
    function obj=get_selInd_table_(obj)
        %Blk key={'mode','lvlInd','blk','trl','intrvl','cmpInd','cmpNum'};
        %obj.lvlIndKey=[ ...
        %              'subInd','lvlInd' ...
        %              repmat({'stdInd'},1,size(stdRC,2)) ...
        %              repmat({'cmpInd'},1,size(cmpRC,2)) ...
        %            ];...
        b=obj.blkOpts;
        repeats=b.repeats;
        mirror=b.mirror;
        binVals=obj.get_src_column('B');
        bins=[obj.ptchsOpts.bins{:}];
        binVals=vertcat(binVals{:});

        nLvl=prod(obj.nLvlPerDim);
        nModes=numel(b.modes);
        nBlk=nLvl.*b.nBlkPerLvl; %50
        nTrlPerBlk=b.nTrlPerLvl./b.nBlkPerLvl; % 180
        nTrlPerMode=nTrlPerBlk*nBlk;
        nTrl=nTrlPerMode*nModes;
        nIntrvlPerMode=nTrlPerMode*b.nIntrvlPerTrl;
        nIntrvlAll=nTrl*b.nIntrvlPerTrl;
        nIntrvlPerBin=nIntrvlAll/numel(bins);


        for i = 1:length(bins)
            obj=bin_table_fun(obj,bins(i),binVals,repeats,mirror,nIntrvlPerBin);
        end

        function obj=bin_table_fun(obj,bin,binVals,srepeats,mirror,nIntrvlPerBin)

            % TODO
            % bins
            % and/or?

            binInd=find(ismember(binVals,bin));
            nStm=numel(binInd);

            bReplace=false;
            if nStm < nIntrvlPerBin
                bReplace=true;
                disp([ 'Setting replace to true. N aviablable stim in bin ' num2str(bin) ' = ' num2str(nStm) '. Req = ' num2str(nIntrvlPerBin) ]);
            end

            if isempty(repeats)
                nUnq=1;
                uInds=ones(nIntrvlPerBin,1);
            else
                col=obj.get_block_column(repeats);
                [~,~,uInds]=unique(col(binInd),'rows');
                nUnq=(max(uInds));
            end


            obj.selInd=zeros(nIntrvlAll,1);
            rng(b.sd);
            for i = 1:nUnq
                ind=(uInds==i);
                K=sum(ind);
                obj.selInd(ind)=datasample(binInd,K,'Replace',bReplace);
            end

            if isempty(mirror)
                return
            else
                col=obj.get_block_column(mirror);
                [~,~,uInds]=unique(col,'rows');
                nUnq=(max(uInds));
            end

            copyInd=(uInds==1);
            for i = 2:nUnq
                ind=(uInds==i);
                K=sum(ind);
                obj.selInd(ind)=obj.selInd(copyInd);
            end
        end
    end

end
methods(Static,Access=private)
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
end
methods(Static)
    function [table,key]=get_blk_table(nDim,nLvlPerDim,nCmpPerLvl,modes,nBlkPerLvl,nTrlPerLvl,nIntrvlPerTrl,sd)

        nLvl=prod(nLvlPerDim);
        nCmp=prod(nCmpPerLvl);
        nModes=numel(modes);
        nBlk=nModes.*nLvl.*nBlkPerLvl; %50
        nTrlPerBlk=nTrlPerLvl./nBlkPerLvl; % 180

        table=distribute(modes,1:nLvl,1:nBlkPerLvl,1:nTrlPerBlk,1:nIntrvlPerTrl); % 18000

        rng(sd);

        cmpNum=Blk.get_cmp_num(nTrlPerBlk,nBlk,nIntrvlPerTrl);
        cmpInd=Blk.get_cmpInd(nCmpPerLvl,nTrlPerBlk,nBlk,nIntrvlPerTrl);

        table=[table cmpInd cmpNum];
        key={'mode','lvlInd','blk','trl','intrvl','cmpInd','cmpNum'};
    end
    function c=get_cmpInd(nCmpPerLvl,nTrlPerBlk,nBlk,nIntrvlPerTrl)
        % which comparison to use
        nIntrvlAll=nTrlPerBlk*nBlk*nIntrvlPerTrl;
        c= zeros(nIntrvlAll,length(nCmpPerLvl));
        for i = 1:length(nCmpPerLvl)
            c(:,i)=cmp_fun(nCmpPerLvl(i),nTrlPerBlk,nBlk,nIntrvlPerTrl);
        end
        function c=cmp_fun(nCmpPerLvl,nTrlPerBlk,nBlk,nIntrvlPerTrl)
            c=repmat(repelem(1:nCmpPerLvl,1,nTrlPerBlk/nCmpPerLvl),nBlk,1);
            c=transpose(Blk.shuffle_within_rows_(c));
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
        c=transpose(Blk.shuffle_within_rows_(c));
        c=c(:);
    end
    function out=parse_blkOpts(Opts)
        P=Blk.get_blk_parseOpts();
        out=parse([],Opts,P);
    end
    function opts=parse_ptchsOpts(opts)
        % XXX TODO fill in missing params if not exist
        %opts=parse([],Opts,P);

        P=Blk.get_ptchs_parseOpts();
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
            elseif endsWith(fld,'XYdeg')
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

    function [table,key,nRow,nCol,nAisle, stdRC,cmpRC,cndRC]=get_opts_tables(opts)
        P=Blk.get_ptchs_parseOpts();
        flds=P(:,1);

        [nRow,nCol,nAisle]=get_dimensions_fun(opts,flds);
        [stdRC,cmpRC,cndRC]=get_RC(nCol,nAisle);
        [table,key]=expand_fun(opts,flds,stdRC,cmpRC);

        function [table,key]=expand_fun(opts,flds,stdRC,cmpRC,cndRC)
            flds(ismember(flds,[{'linked','dims'} opts.dims]))=[];
            N=size(stdRC,1);
            M=length(flds);
            ndim=size(stdRC,2);

            table=cell(N, M);
            dimTable=cell(N,ndim*2);
            for j = 1:M
                fld=flds{j};
                if ismember(fld,horzcat(opts.linked{:}))
                    table(:,j)=expand_cell_fun(opts,fld);
                elseif ismember(fld,opts.dims)
                    ind=find(ismember(opts.dims,fld));
                    [dimTable(:,ind), dimTable(:,ind+ndim)] =...
                        expand_dim_fun(opts,fld,ind,stdRC,cmpRC,N);
                elseif iscell(opts.(fld))
                    error('Cells are reserved for linked params');
                else
                    table(:,j)=expand_single_fun(opts,fld,N);
                end
            end
            table=[dimTable table];
            f=transpose(flds);
            f(ismember(f,opts.dims))=[];
            key=[ repmat({'std'},1,ndim) repmat({'cmp'},1,ndim) f];
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

        function [stdRC,cmpRC,cndRC]=get_RC(nCol,nAisle)
        %std{:} cmp{:}
            C=cellfun(@(x) 1:x, num2cell(nCol),UO,false);
            A=cellfun(@(x) 1:x, num2cell(nAisle),UO,false);
            RC=distribute(distribute(C{:}),distribute(A{:}));
            h=size(RC,2)/2;
            stdRC=RC(:,1:h);
            cmpRC=RC(:,h+1:end);
            [~,~,u]=unique(stdRC,'rows');
            cndRC=[transpose(1:length(u) ) u];
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
        P=Blk.get_blk_parseOpts();
        flds=P(:,1);
    end
    function flds=get_ptchs_flds();
        P=Blk.get_ptchs_parseOpts();
        flds=P(:,1);
    end
    function P=get_blk_parseOpts()
        P={...
              'dims',[],[] ...
              ;'repeats',[],'' ...
              ;'mirror',[],'' ...
              ;'modes',[],'isallint' ...
              ;'nBlkPerLvl',[],'isint' ...
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

     
          };
    end
end
end
