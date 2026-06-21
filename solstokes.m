% Compute Stokes solution for the lid driven cavity problem
copyfile('./stokes_flow/test_problems/tightcavity_flow.m', './stokes_flow/specific_flow.m');
copyfile('./stokes_flow/test_problems/zero_bc.m', './stokes_flow/stream_bc.m');
square_stokes;