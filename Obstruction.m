classdef Obstruction
    properties
        origin (1,2)
        x_size (1,1)
        y_size (1,1)
        BC = 0
        diag_up
        diag_dn
    end
    methods
        function obj = Obstruction(o, x, y, bc)
            obj.origin = o;
            obj.x_size = x;
            obj.y_size = y;
            obj.BC = bc;
            m = y/x;
            obj.diag_up = @(x) m.*(x - o(1)) + o(2);
            obj.diag_dn = @(x) -m.*(x - o(1)) + o(2) + y;
        end
        function [IN, INL, INR, INU, IND] = isInside(p, obj)
            IN = p(1,:) >= obj.origin(1) & p(1,:) <= (obj.origin(1) + obj.x_size) ...
                & p(2,:) >= obj.origin(2) & p(2,:) <= (obj.origin(2) + obj.y_size);
            INL = IN & p(2,:) >= obj.diag_up(p(1,:)) & p(2,:) <= obj.diag_dn(p(1,:));
            INR = IN & p(2,:) <= obj.diag_up(p(1,:)) & p(2,:) >= obj.diag_dn(p(1,:));
            INU = IN & p(2,:) >= obj.diag_up(p(1,:)) & p(2,:) >= obj.diag_dn(p(1,:));
            IND = IN & p(2,:) <= obj.diag_up(p(1,:)) & p(2,:) <= obj.diag_dn(p(1,:));
        end
    end
end