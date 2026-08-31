classdef I18n < handle
% I18n  Bilingual string table loaded from zh.json / en.json.
    properties (Access = private)
        tables (1,1) struct = struct('zh', struct(), 'en', struct())   % .zh, .en  each a struct(key->text)
    end
    properties
        lang (1,:) char = 'zh'
    end
    methods
        function obj = I18n(resourceDir)
            arguments
                resourceDir (1,:) char
            end
            obj.tables.zh = obj.readJson(fullfile(resourceDir, 'zh.json'));
            obj.tables.en = obj.readJson(fullfile(resourceDir, 'en.json'));
        end
        function setLang(obj, lang)
            arguments
                obj
                lang (1,:) char {mustBeMember(lang, {'zh','en'})}
            end
            obj.lang = lang;
        end
        function s = t(obj, key)
            arguments
                obj
                key (1,:) char
            end
            tbl = obj.tables.(obj.lang);
            if isfield(tbl, key)
                s = tbl.(key);
            else
                s = key;   % fallback: show the key itself
            end
        end
    end
    methods (Static, Access = private)
        function tbl = readJson(path)
            % Read explicitly as UTF-8 so Chinese text is not corrupted by
            % the system default charset (e.g. GBK on Windows) before decode.
            fid = fopen(path, 'r', 'n', 'UTF-8');
            if fid < 0
                error('I18n:NoFile', 'i18n file not found: %s', path);
            end
            txt = fread(fid, inf, '*char')';
            fclose(fid);
            tbl = jsondecode(txt);
        end
    end
end
