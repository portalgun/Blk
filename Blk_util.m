classdef Blk_util < handle
methods
end
methods(Static)
    function [t,dims] = make_lookup_tables(BC)
        lookup=Blk.make_lookup_struct(BC);
        [t,dims]=Blk.lookup_to_tables(lookup);
    end
    function lookup = make_lookup_struct(obj)
        lookup=struct;
        if isa(obj,'Blk_con')
            lookup.dims=obj.dims;
            lookup.cnd=obj.cndLookup;
            lookup.cndKey=obj.cndKey;
            lookup.lvl=obj.lvlLookup;
            lookup.lvlKey=obj.lvlKey;
            lookup.cmp=obj.cmpLookup;
            lookup.cmpKey=obj.cmpKey;
        else
            lookup.dims=obj.dims;
            lookup.cnd=obj.lookup.cnd.TABLE;
            lookup.cndKey=obj.lookup.cnd.KEY;
            lookup.lvl=obj.lookup.lvl.TABLE;
            lookup.lvlKey=obj.lookup.lvl.KEY;
            lookup.cmp=obj.lookup.cmp.TABLE;
            lookup.cmpKey=obj.lookup.cmp.KEY;
        end
    end
    function [new,dims]=lookup_to_tables(lookup)
        new=struct();
        new.cnd=Table(lookup.cnd,lookup.cndKey);
        new.lvl=Table(lookup.lvl,lookup.lvlKey);
        new.cmp=Table(lookup.cmp,lookup.cmpKey);
        dims=lookup.dims;
    end
    % LOAD
    function [t,newInd]=load_blk(alias)
        dire=Blk.getDir(alias);
        fname=[dire 'blk.mat'];
        S=load(fname);
        types=repmat({'double'}, 1,length(S.key));
        t=Table(S.table,S.key,types);
        if isfield(S,'newInd')
            newInd=S.newInd;
        else
            newInd=[];
        end
    end
    function t=load_opts(alias)
        S=Blk.load_opts_raw(alias);
        t=Table(S.table,S.key);
    end
    function [t,dims]=load_lookup(alias)
        dire=Blk.getDir(alias);
        S=load([dire 'lookup.mat']);
        dims=S.lookup.dims;
        t=Blk.lookup_to_tables(S.lookup);
    end
    function opts=load_opts_raw(alias)
        dire=Blk.getDir(alias);
        opts=load([dire 'opts.mat']);
    end
    % DIR
    function dire=getDir(alias)
        rt=Env.var('root');
        if isempty(rt)
            rt=Env.var('BLK__ROOT');
        end
        dire=[rt alias filesep];
    end
    % FNAME
    function fname=get_opts_table_fname(alias)
        fname=[Blk.getDir(alias) 'opts'];
    end
    function fname=get_blk_table_fname(alias)
        fname=[Blk.getDir(alias) 'blk'];
    end
    function fname=get_lvlInd_table_fname(alias)
        fname=[Blk.getDir(alias) 'table'];
    end
end
end
