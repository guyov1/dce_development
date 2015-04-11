function Filter = Adjusted_Larsson_Filter_no_Ve(time_vec_min, F, Vb, E)

% No extra-vascular - Uptake Model
Tp     = (Vb/F) * (1-E);
Ktrans = E*F;
IRF    = F * exp(-(1/Tp) * time_vec_min) + Ktrans*(1-exp(-(1/Tp) * time_vec_min));
Filter = IRF;

end