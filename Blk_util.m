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
%% CONVERSIONS

    %% CND
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
