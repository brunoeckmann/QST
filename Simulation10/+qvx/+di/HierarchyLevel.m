classdef HierarchyLevel
    properties
        hierarchyType = '';
        level = 0;
        addNonNegativity = 0;
    end
    methods
        function HL = HierarchyLevel(hierarchyType, level, addNonNegativity)
        % do not use the constructor, but the static methods
        % HierarchyLevel.NPA(n) or HierarchyLevel.local(n)
            HL.hierarchyType = hierarchyType;
            HL.level = level;
            HL.addNonNegativity = addNonNegativity;
        end
    end
    methods (Static)
        function HL = NPA(level)
        % Returns the NPA hierarchy level as defined in Navascues, Pironio, Acin, 2008
        % For the level 1, does not include the positivity constraints
            HL = qvx.di.HierarchyLevel('NPA', level, false);
        end
        function HL = NPAplus(level)
        % Similar to HierarchyLevel.NPA, but adds nonnegativity constraints
        % when level = 1
            HL = qvx.di.HierarchyLevel('NPA', level, level == 1);
        end
        function HL = local(level)
        % Returns the local level variant described in Moroder 2013
            HL = qvx.di.HierarchyLevel('local', level, false);
        end
        function HL = almostQuantum
            HL = qvx.di.HierarchyLevel.local(1);
        end
    end
end
