classdef spline1
    %1SPLINE a simple class for linear interpolants. Basic functionality.
    %No documentation is provided. Written by K. Church, 2021.
    
    properties (Access=protected)
        a           % degree 1 coefficient vector(s)
        b           % degree 0 coefficient vector(s)
        domain      % closed interval(s)
        pieces      % number of pieces
        segments    % number of segments, per piece (assumed equal)
    end
    
    methods
        function s = spline1(a,b,domain)
            %1SPLINE Construct a linear spline from its coefficients on
            %each sub-interval. 
            %Input: a ------ degree 1 coefficient vector
            %       b ------ degree 0 coefficient vector
            %       domain - domain of the spline
            if size(a,1)~=size(b,1)
                error('Coefficients must have the same number of rows.');
            elseif size(a,2)~=size(b,2)
                error('Coefficients must have the same number of columns.');
            elseif size(domain,2)~=size(b,2)
                error('Number of domains does not match number of columns of coefficients');
            elseif size(domain,1)~=2
                error('Domain must be 2xN matrix.');
            end
            chkdom = 1;
            for k=1:size(a,2)-1
                chkdom = min(chkdom,domain(2,k)==domain(1,k+1));
            end
            if chkdom == 0
                error('Incompatible domains');
            end
            s.a = a;
            s.b = b;
            s.domain = domain;
            s.segments = size(a,1);
            s.pieces = size(a,2);
        end
        
        function disp(s)
            if s.pieces>1
                fprintf('Linear spine; %s pieces each with %s sub-intervals on [%s,%s]. Coefficients suppressed. \n',...
                    num2str(s.pieces),num2str(s.segments),...
                    num2str(s.domain(1,1)),num2str(s.domain(end,end)));
            else
                fprintf('Linear spline on %s sub-intervals of [%s,%s].\n',...
                num2str(s.segments),num2str(s.domain(1)),...
                num2str(s.domain(2)));
                fprintf('[a|b] = \n');
                disp([num2str(s.a),repmat(' | ',[s.segments,1]),num2str(s.b)]);
            end
        end
        
        function plot(s)
            M = s.segments; dom = s.domain;
            for k=1:s.pieces
                x = dom(1,k) + (dom(2,k)-dom(1,k))*(0:M).'/M;
                y = eval1spline(get_subpiece(s,k),x);
                plot(x,y,'k');
                hold on
            end
            hold off
        end
        
        function a = get_a(s)
            a = s.a;
        end
        
        function b = get_b(s)
            b = s.b;
        end
        
        function seg = get_segments(s)
            seg = s.segments;
        end
        
        function dom = get_dom(s)
            dom = s.domain;
        end
        
        function pieces = get_pieces(s)
            pieces = s.pieces;
        end
        
        function s_k = get_subpiece(s,k)
           s_k = spline1(s.a(:,k),s.b(:,k),s.domain(:,k));
        end
        
        function [s,w] = fuse(s1,s2,flag_it)
           s = spline1([s1.a,s2.a],[s1.b,s2.b],[s1.domain,s2.domain]);
           M = get_segments(s1);
           s1_end = s1.b(end,s1.pieces)+s1.a(end,s1.pieces)/M;
           s2_start = s2.b(1,1);
           if nargin<3
               flag_it = 1;
           end
           w = abs(s1_end - s2_start);
           if abs(s1_end - s2_start)>1E-14 && flag_it == 1
               warning('gap between segments to be fused is greater than 1E-14; possible discontinuity.');
           end
        end
        
        function s = shift(s,dom_shift)
            s.domain = s.domain + dom_shift;
        end
    end
end