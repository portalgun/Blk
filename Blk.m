classdef Blk < handle & Blk_util & Blk_verify
properties
    alias

    blk
    opts
end
properties(Hidden)
    dims
    lookup
    lims={}
end
methods
    function obj=Blk(defName,replaceBinds,bSave)
        if nargin < 1
            return
        end
        if ~exist('replaceBinds','var')
            replaceBinds=[];
        end
        if ~exist('bSave','var') || isempty(bSave)
            bSave=1;
        end
        BC=Blk_con(defName,replaceBinds,bSave);
        obj.alias=BC.alias;

        obj.blk=Table(BC.blkTable, BC.blkKey);
        obj.opts=struct();
        obj.opts.table=BC.optsTable;
        obj.opts.key=BC.optsKey;
        %obj.opts=Table(BC.optsTable, BC.optsKey); % XXX
        [obj.lookup,obj.dims]=Blk.make_lookup_tables(BC);
    end
    function new=select_block(obj,mode,lvlInd,blkNum)
        new=copyObj(obj);
        new.dims=obj.dims;
        new.lookup=obj.lookup;
        new.blk=obj.blk('mode',mode,'blk',blkNum,'lvlInd',lvlInd);
        inds=new.blk.unique('cndInd');
        new.opts.table=obj.opts.table(inds,:);
        new.lims={'mode', mode ;'blk', blkNum; 'lvlInd', lvlInd};
        %new.opts.table
    end
    function S=ret_opts_struct(obj)
        % EG ptchs.ptchOpts
        S=struct;
        n=length(obj.opts.key);
        bExpand=false(n,1);
        bExpandD=false(n,1);
        bExpandL=false(n,1);
        bExpandDL=false(n,1);
        for i = 1:n
            fld=obj.opts.key{i};
            col=obj.get_opts_column(fld);

            if isuniform(col)
                S.(fld)=col{1};
                if iscell(S.(fld)) && numel(S.(fld)) ==1
                    S.(fld)=S.(fld){1};
                elseif ismember(fld,obj.dims)
                    bExpandD(i)=true;
                elseif iscell(S.(fld)) && all(cellfun(@isnumeric,S.(fld)))
                    bExpand(i)=true;
                end
            else
                S.(fld)=col;
                if ismember(fld,obj.dims)
                    bExpandDL(i)=true;
                elseif iscell(S.(fld){1})
                    bExpandL(i)=true;
                end
            end

        end
        flds=obj.opts.key(bExpand);
        fldsD=obj.opts.key(bExpandD);
        fldsL=obj.opts.key(bExpandL);
        fldsDL=obj.opts.key(bExpandDL);
        if ~isempty(fldsD)
            S=obj.expand_fun_dims(S,fldsD);
        end
        if ~isempty(fldsDL)
            S=obj.expand_fun_dims_long(S,fldsDL);
        end

    end
%% COL
    function col=get_opts_column(obj,key)
        ind=ismember(obj.opts.key,key);
        col=obj.opts.table(:,ind);
    end
    function col=get_block_column(obj,key)
        if ~iscell(key)
            key={key};
        end
        blkKey=obj.blk.get_key();
        gd=ismember(key,blkKey);
        if all(gd)
            col=obj.blk(key{:}).ret();
            return
        end

        col=zeros(size(obj.blk,1), numel(key));
        ind=ismember(blkKey,key(gd));
        if any(gd)
            col(:,gd)=obj.blk(key{ind}).ret();
        end
        col(:,~gd)=obj.get_block_column_dim_lvl(key(~gd));
    end
    function col=get_block_column_dim_lvl(obj,name)
        col=zeros(size(obj.blk,1), numel(name));
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
            a=obj.blk.lvlInd_to_lvlRC(lvls);
            col(:,i)=a(:,d);
        end
    end
%% CONVERSIONS
    function RC=lvlInd_to_lvlRC(obj,ind)
        RC=obj.lookup.lvl.ret();
        RC=RC(ind,2:3);
    end
    function RC=lvlRC_to_lvlX(obj,RC)
    end
    function lvl=cndInd_to_lvlInd(obj,cndInd)
        lvl=obj.lookup.cndr(cndInd,'lvlInd').ret();
    end
    function cmps=cndInd_to_cmpInd(obj,cndInd)
        cmps=obj.lookup.cnd(cndInd,'cmpInd').ret();
    end

    function cnd=lvlInd_cmpInd_to_cndInd(obj,lvl,cmp)
        nStd=size(obj.lookup.lvl,1);
        nCmp=size(obj.lookup.cmp,1);

        cnd=sub2ind([nCmp nStd],cmp,lvl);
    end
    function stdX=cndInd_to_stdX(obj,cndInd)
        col=obj.get_opts_column(obj.dims);
        if ~iscell(cndInd)
            cndInd=num2cell(cndInd);
        end
        t=obj.blk('cndInd').unique();
        cndInd=cellfun(@(x) find(t==x,1,'first'), cndInd);
        %cndInd=obj.blk.find('cndInd',cndInd);
        stdX=zeros(numel(cndInd),numel(obj.dims));
        for i = 1:size(col,2)
            % XXX
            s=vertcat(col{cndInd,i});
            stdX(:,i)=cell2mat(s(:,1));
        end
    end
    function cmpX=cndInd_to_cmpX(obj,cndInd,cmpNum)
        if ~exist('cmpNum','var') || isempty(cmpNum)
            cmpNum=1;
        end
        if ~iscell(cndInd)
            cndInd=num2cell(cndInd);
        end
        t=obj.blk('cndInd').unique();
        cndInd=cellfun(@(x) find(t==x,1,'first'), cndInd);

        ci=cmpNum+1;
        col=obj.get_opts_column(obj.dims);
        cmpX=zeros(numel(cndInd),numel(obj.dims));
        for i = 1:size(col,2)
            s=vertcat(col{cndInd,i});
            cmpX(:,i)=cell2mat(s(:,ci));
        end
    end
    function cndInd=trial_to_cndInd(obj,trials)
        if ~iscell(trials)
            trials=num2cell(trials);
        end
        cndInd=obj.blk('trl',trials{:},'intrvl',1,'cndInd').ret();
    end
    function stdX=trial_to_stdX(obj,trial)
        cndInd=obj.trial_to_cndInd(trial);
        stdX=obj.cndInd_to_stdX(cndInd);
    end
    function cmpX=trial_to_cmpX(obj,trial,cmpNum)
        if ~exist('cmpNum','var')
            cmpNum=[];
        end
        cndInd=obj.trial_to_cndInd(trial);
        cmpX=obj.cndInd_to_cmpX(cndInd,cmpNum);
    end
    function stmInd=trial_to_stmInd(obj,trls)
        if ~iscell(trls);
            trls=num2cell(trls);
        end
        stmInd=obj.blk.find('trl',trls{:});
    end
    function [trl,intrvl]=stmInd_to_trial_interval(obj,stmInd)
        if ~iscell(stmInd);
            stmInd=num2cell(stmInd);
        end
        out=obj.blk(stmInd{:},'trl','intrvl').ret();
        trl=out(:,1);
        intrvl=out(:,2);
    end
    function stmInd=trial_intrvl_to_stmInd(obj,trls,intrvls)
        if ~iscell(trls);
            trls=num2cell(trls);
        end
        if ~iscell(intrvls);
            intrvls=num2cell(intrvls);
        end
        stmInd=obj.blk.find('trl',trls{:},'intrvl',intrvls{:});
    end

%% GET
    function bMotion=get_bMotion(obj)
        bMotion=all(~cellfun(@isempty,obj.get_opts_column('speed')));
    end
    function nTrl=get_nTrial(obj)
        nTrl=numel(obj.blk.unique('trl'));
    end
    function nStm=get_nStm(obj)
        nStm=sum(obj.blk('P').ret() > 0);
    end
    function nIntrvl=get_nIntrvl(obj)
        nIntrvl=numel(obj.blk.unique('intrvl'));
    end
    function nCmp=get_nCmp(obj)
        nCmp=sum(obj.blk.unique('cmpNum')>0);
    end
    function cmpInt=get_cmpIntrvl(obj,trl,cmpNum)
        if ~exist('cmpNum') || isempty(cmpNum)
            cmpNum=1;
        end
        if ~exist('trl','var') || isempty(trl)
            nTrl=1:obj.get_nTrial();
            trl=num2cell(nTrl);
        end
        if ~iscell(trl)
            trl=num2cell(trl);
        end
        cmpInt=obj.blk('trl',trl{:},'cmpNum',cmpNum,'intrvl').ret();
    end

    function stdX=get_stdX(obj,trl)
        if ~exist('trl','var') || isempty(trl)
            nTrl=obj.get_nTrial();
            trl=1:nTrl;
        end
        stdX=obj.trial_to_stdX(trl);
    end
    function cmpX=get_cmpX(obj,trl,cmpNum)
        if ~exist('cmpNum','var')
            cmpNum=[];
        end
        if ~exist('trl','var') || isempty(trl)
            nTrl=obj.get_nTrial();
            trl=1:nTrl;
        end

        cmpX=obj.trial_to_cmpX(trl,cmpNum);
    end
    function stdXunq=get_stdX_unq(obj)
        stdX=obj.get_stdX();
        stdXunq=unique(stdX,'rows');
    end
    function cmpXunq=get_cmpX_unq(obj,cmpNum)
        if ~exist('cmpNum','var')
            cmpNum=[];
        end
        cmpX=obj.get_cmpX();
        cmpXunq=unique(cmpX,'rows');
    end
    function answ=get_correct(obj)
        nCmp=obj.get_nCmp();

        if nCmp==1
            answ=obj.get_correct_2IFC();
        end
    end
    function answ=get_correct_2IFC(obj);
        stdX=obj.get_stdX();
        cmpX=obj.get_cmpX([],1);
        cmpIntrvl=obj.get_cmpIntrvl();

        answ=Rsp.get_correct_2IFC([],stdX,cmpX,cmpIntrvl);
    end
%% RET
    function S=ret_blk_struct(obj)
    % EG expTracker
        S=struct();

        %numel(obj.blk.unique('cndInd')) % 125
        S.nStd=numel(obj.blk.unique('lvlInd'));  % 25 | 5
        S.nCmp=numel(obj.blk.unique('cmpInd'));  % 5  | 9
        S.nBlk=numel(obj.blk.unique('blk'));     % 5  | 5
        S.nTrlPerBlk=numel(obj.blk.unique('trl'));     %   | 900
        S.nTrlPerLvl=S.nTrlPerBlk*S.nBlk;
        S.nTrl=S.nTrlPerLvl*S.nCmp;
        S.Xname=obj.dims;

        S.methodVars=obj.get_method_vars();

    end
    function m=get_method_vars(obj)
        m=struct();
        [m.stdXunqAll,m.cmpXunqAll]=obj.get_stds_cmps_unq();
        m.minDist=[];

        %obj.blk

    end
    function [stds,cmps]=get_stds_cmps_unq(obj)
        Cnds=obj.lookup.cnd.ret();
        Lvls=obj.lookup.lvl.ret();
        Cmps=obj.lookup.cmp.ret();

        stds=cell(length(obj.dims),1);
        cmps=cell(length(obj.dims),1);
        inds=zeros(length(obj.dims),1);
        for i = 1:length(obj.dims)
            dim=obj.dims{i};
            nStd=numel(unique(Lvls(:,i+1)));
            nCmp=numel(unique(Cmps(:,i+1)));
            stds{i}=cell(nStd,nCmp);
            cmps{i}=cell(nStd,nCmp);
            inds(i)=find(ismember(obj.opts.key,obj.dims{i}));
        end
        for i = 1:size(Cnds,1)
            lvlInd=Cnds(i,2);
            cmpInd=Cnds(i,3);
            lvlRC=Lvls(lvlInd,2:end);
            cmpRC=Cmps(cmpInd,2:end);

            r=obj.opts.table(i,:);
            for j = 1:length(obj.dims)
                dim=obj.dims{j};
                v=r{inds(j)};

                stds{j}( lvlRC(j), cmpRC(j) ) =v(1);
                cmps{j}( lvlRC(j), cmpRC(j) ) =v(2:end);
            end
        end
        for i = 1:length(obj.dims)
            dim=obj.dims{i};
            if all(cellfun(@isnumeric,stds{i}))
                stds{i}=cell2mat(stds{i});
                if isuniform(stds{i},2)
                    stds{i}=stds{i}(:,1);
                end
            end
            if all(cellfun(@isnumeric,cmps{i}))
                cmps{i}=cell2mat(cmps{i});
            end
        end

    end
    function cmps=get_cmp_unq(obj)
        cmps=obj.lookup.cmp(:,2:end);
    end
end
methods(Access=private)
    function S=expand_fun_dims(obj,S,fldsD)
        % cmpnum
        cmpNum=obj.blk('cmpNum').ret()+1;
        for i = 1:length(fldsD)
            fld=fldsD{i};
            val=transpose(S.(fld)(cmpNum));
            if all(cellfun(@isnumeric,val))
                S.(fld)=cell2mat(val);
            else
                S.(fld)=val;
            end
        end
    end
    function S=expand_fun_dims_long(obj,S,fldsDL)
        % cndInd
        cmpNum=obj.blk('cmpNum').ret()+1;
        %cndInd=obj.blk('cndInd').ret();
        [~,~,cndInd]=obj.blk('cndInd').unique('cndInd'); % XXX?
        sz=[length(cndInd), max(cmpNum)];
        ind=sub2ind(sz, transpose(1:sz(1)), cmpNum);
        for i = 1:length(fldsDL)
            fld=fldsDL{i};

            val=S.(fld)(cndInd);
            val=vertcat(val{:});
            val=val(ind);
            if all(cellfun(@isnumeric,val))
                S.(fld)=cell2mat(val);
            else
                S.(fld)=val;
            end
        end
    end
end
methods(Static)
    function obj=load(alias)
        obj=Blk();
        obj.alias=alias;
        obj.blk=Blk.load_blk(alias);
        obj.opts=Blk.load_opts_raw(alias); % XXX change later
        [obj.lookup,obj.dims]=Blk.load_lookup(alias);
    end
end
end
