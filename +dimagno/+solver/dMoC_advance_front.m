%----------------------------------------------------------------------
%{

This function receives a integration front_2d (Fold) and propagates it
downstream once (Fnew).

%}
%----------------------------------------------------------------------

function [Fnew] = dMoC_advance_front(data,Fold)
  
Fsub = dimagno.front;
Fnew = dimagno.front;

%% Subfront    
for i = 1:Fold.n_points-1         
    p = dimagno.solver.dMoC_inner(data,i+1,i,Fold);
    Fsub.setpoint(p,i); clear p;   
end

%% Calculate inner
for i = 1:Fsub.n_points-1 
    p = dimagno.solver.dMoC_inner(data,i+1,i,Fsub);
    Fnew.setpoint(p,i); clear p;   
end

%% Calculate axis       
Fsub.setpoint(Fold.point(1),'prepend');
p = dimagno.solver.dMoC_axis(data,2,1,Fsub);
Fnew.setpoint(p,'prepend'); clear p;       

%% Calculate border point    
Fsub.setpoint(Fold.point(Fold.n_points),'append');
p = dimagno.solver.dMoC_border(data,Fsub.n_points,Fsub.n_points-1,Fsub);
Fnew.setpoint(p,'append'); clear p;
     
