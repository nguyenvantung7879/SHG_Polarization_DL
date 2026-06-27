function [Fiber_Angle_Predict, Peptide_Pitch_Angle_Predict, Time_predict] = predict_image(net, I_SHG, Threshold)

   
    % Optional: threshold mask step if needed
    if canUseGPU
        I_SHG = gpuArray(single(I_SHG));
    else
        I_SHG = single(I_SHG);
    end
    
    tic;
    % === Predict angle map from SHG stack ===
    Angle_Predict = activations(net, I_SHG, 'Angle_Predict', ...
                                 'ExecutionEnvironment', 'gpu');

    Time_predict = toc;
    % === Extract angle channels and apply correction ===
    Fiber_Angle_Predict = Angle_Predict(:, :, 1);
    Peptide_Pitch_Angle_Predict = Angle_Predict(:, :, 2);


end
