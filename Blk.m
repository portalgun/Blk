classdef Blk < handle & Blk_util & Blk_verify
properties
    alias

    blk
    opts
end
properties(Hidden)
    dims
    lookup
    lims={};
end
methods
    function obj=Blk(defName)
        if nargin < 1
            return
        end
        BC=Blk_con(defName);
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
    function col=get_opts_column(obj,key)
        ind=ismember(obj.opts.key,key);
        col=obj.opts.table(:,ind);
    end
    function S=ret_blk_struct(obj)
    % EG expTracker
        S=struct();

        %numel(obj.blk.unique('cndInd')) % 125
        S.nStd=numel(obj.blk.unique('lvlInd'));  % 25 | 5
        S.nCmp=numel(obj.blk.unique('cmpInd'));  % 5  | 9
        S.nBlk=numel(obj.blk.unique('blk'));     % 5  | 5
        S.nTrl=[];
        S.nTrlPerBlk=numel(obj.blk.unique('trl'));     %   | 900
        S.nTrl=S.nTrlPerBlk*S.nBlk;
        S.nTrlPerLvl=S.nTrl/S.nCmp; % cmp is lvl?
        S.Xname=obj.dims;

        S.methodVars=obj.get_method_vars();

    end
    function m=get_method_vars(obj)
        m=struct();
        [m.stdXunqAll,m.cmpXunqAll]=get_stds_cmps_unq();
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
        cndInd=obj.blk('cndInd').ret();
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
