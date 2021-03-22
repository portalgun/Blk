classdef Blk_util < handle
methods
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
%% UTIL
    % XXX MOVE TO CON
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
    function col=get_src_column(obj,key)
        ind=ismember(obj.srcKey,key);
        col=obj.srcTable(:,ind);
    end
%% CONVERSIONS

    %% CND
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
        RC=[obj.lvlLookup(ind,2:3)];
    end
    function RC=cmpInd_to_cmpRC(obj,ind)
        RC=[obj.cmpLookup(ind,2:3)];
    end
end
methods(Static)
    function [t,dims] = make_lookup_tables(BC)
        lookup=Blk.make_lookup_struct(BC);
        [t,dims]=Blk.lookup_to_tables(lookup);
    end
    function lookup = make_lookup_struct(BC)
        lookup=struct;
        lookup.dims=BC.dims;
        lookup.cnd=BC.cndLookup;
        lookup.cndKey=BC.cndKey;
        lookup.lvl=BC.lvlLookup;
        lookup.lvlKey=BC.lvlKey;
        lookup.cmp=BC.cmpLookup;
        lookup.cmpKey=BC.cmpKey;
    end
    function [new,dims]=lookup_to_tables(lookup)
        new=struct();
        new.cnd=Table(lookup.cnd,lookup.cndKey);
        new.lvl=Table(lookup.lvl,lookup.lvlKey);
        new.cmp=Table(lookup.cmp,lookup.cmpKey);
        dims=lookup.dims;
    end
    % LOAD
    function t=load_blk(alias)
        dire=Blk.get_dir(alias);
        fname=[dire 'blk.mat'];
        S=load(fname);
        t=Table(S.table,S.key);
    end
    function t=load_opts(alias)
        S=Blk.load_opts_raw(alias);
        t=Table(S.table,S.key);
    end
    function [t,dims]=load_lookup(alias)
        dire=Blk.get_dir(alias);
        S=load([dire 'lookup.mat']);
        dims=S.lookup.dims;
        t=Blk.lookup_to_tables(S.lookup);
    end
    function opts=load_opts_raw(alias)
        dire=Blk.get_dir(alias);
        opts=load([dire 'opts.mat']);
    end
    % DIR
    function dire=get_dir(alias)
        dire=[dbDirs('blk') alias filesep];
    end
    % FNAME
    function fname=get_opts_table_fname(alias)
        fname=[Blk.get_dire(alias) 'opts'];
    end
    function fname=get_blk_table_fname(alias)
        fname=[Blk.get_dire(alias) 'blk'];
    end
    function fname=get_lvlInd_table_fname(alias)
        fname=[Blk.get_dire(alias) 'table'];
    end
end
end
