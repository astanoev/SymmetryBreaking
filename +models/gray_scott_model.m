classdef gray_scott_model < handle
    properties
        state
        D
        par
        F
    end
    
    methods
        function obj = gray_scott_model(state)
            [Du, Dv, F, k] = deal(0.16, 0.08, 0.035, 0.065); % Bacteria 1
            [Du, Dv, F, k] = deal(0.14, 0.06, 0.035, 0.065); % Bacteria 2
            %[Du, Dv, F, k] = deal(0.16, 0.08, 0.060, 0.062); % # Coral
            %[Du, Dv, F, k] = deal(0.19, 0.05, 0.060, 0.062); % # Fingerprint
            %[Du, Dv, F, k] = deal(0.10, 0.10, 0.018, 0.050); % # Spirals
            %[Du, Dv, F, k] = deal(0.12, 0.08, 0.020, 0.050); % # Spirals Dense
            %[Du, Dv, F, k] = deal(0.10, 0.16, 0.020, 0.050); % # Spirals Fast
            %[Du, Dv, F, k] = deal(0.16, 0.08, 0.020, 0.055); % # Unstable
            %[Du, Dv, F, k] = deal(0.16, 0.08, 0.050, 0.065); % # Worms 1
            %[Du, Dv, F, k] = deal(0.16, 0.08, 0.054, 0.063); % # Worms 2
            %[Du, Dv, F, k] = deal(0.16, 0.08, 0.035, 0.060); % # Zebrafish

            %[Du, Dv, F, k] = deal(0.1, 0.08, 0.02, 0.05); % # Spiral waves (pseudo-periodic behavior) 
            %[Du, Dv, F, k] = deal(0.18, 0.09, 0.02, 0.056); % # Chaos (Unstable spot solitons) 
            %[Du, Dv, F, k] = deal(0.01, 0.16, 0.02, 0.05); % # Spiral waves (pseudo-periodic) 
            %[Du, Dv, F, k] = deal(0.18, 0.13, 0.025, 0.056); % # Self-limiting population (borderline chaotic behavior) 
            %[Du, Dv, F, k] = deal(0.14, 0.06, 0.035, 0.065); % # Multiplying spots (static spike solitons) 
            %[Du, Dv, F, k] = deal(0.16, 0.08, 0.035, 0.060); % # Zebrafish (nearly-stable labyrinthine forms resolving to parallel bands) 
            %[Du, Dv, F, k] = deal(0.16, 0.08, 0.05, 0.065); % # Worms (sparse, unbounded stable forms) 
            %[Du, Dv, F, k] = deal(0.19, 0.09, 0.062, 0.062); % # Coral (stable labyrinthine forms) 
            %[Du, Dv, F, k] = deal(0.16, 0.08, 0.03, 0.058); % # Gradient
            obj.state = state;
            obj.D = [Du, Dv];
            obj.par = struct('F',F,'k',k);
        end
        
        function F = reaction(obj, state_inner)
            u = state_inner(1,:,:);
            v = state_inner(2,:,:);
            uvv = u.*v.*v;
            fuv = - uvv +  obj.par.F.*(1-u);
            guv = uvv - (obj.par.F+obj.par.k).*v;
            F = zeros(size(state_inner));
            F(1,:,:) = fuv;
            F(2,:,:) = guv;
        end
    end
end

