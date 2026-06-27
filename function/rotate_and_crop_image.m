function img_rot_scale = rotate_and_crop_image(img, angle_deg)
    % Get original size
    [h, w, ~] = size(img);

    % Rotate image without losing content
    img_rot = imrotate(img, angle_deg, "nearest", "crop");

    % Get size of largest axis-aligned rect within the rotated box
    [maxWidth, maxHeight] = largest_rotated_rect(w, h, angle_deg);

    % Crop around center of the rotated image
    img_rot_crop = crop_around_center(img_rot, maxWidth, maxHeight);

    % Resize back to original size
    img_rot_scale = imresize(img_rot_crop, [256, 256], "nearest");
end

function [maxWidth, maxHeight] = largest_rotated_rect(w, h, angle_deg)

    if angle_deg<0
        angle_deg = 360 - abs(angle_deg);
    end
    % Convert angle from degrees to radians
    angle = deg2rad(angle_deg);

    % Normalize the angle to quadrant
    quadrant = bitand(floor(angle / (pi/2)), 3);
    if bitand(quadrant, 1) == 0
        sign_alpha = angle;
    else
        sign_alpha = pi - angle;
    end
    alpha = mod(mod(sign_alpha, pi) + pi, pi);

    % Compute bounding box width and height
    bb_w = w * cos(alpha) + h * sin(alpha);
    bb_h = w * sin(alpha) + h * cos(alpha);

    % Compute angle gamma based on width and height relationship
    if w < h
        gamma = atan2(bb_h, bb_w); 
    else
        gamma = atan2(bb_w, bb_h); 
    end

    delta = pi - alpha - gamma;

    if w < h
        length = h;
    else
        length = w;
    end

    d = length * cos(alpha);
    a = d * sin(alpha) / sin(delta);

    y = a * cos(gamma);
    x = y * tan(gamma);

    maxWidth = bb_w - 2 * x;
    maxWidth = round(maxWidth);
    maxHeight = bb_h - 2 * y;
    maxHeight = round(maxHeight);
end

function cropped = crop_around_center(image, width, height)
    [imgHeight, imgWidth, ~] = size(image);
    centerX = floor(imgWidth / 2);
    centerY = floor(imgHeight / 2);

    width = min(width, imgWidth);
    height = min(height, imgHeight);

    x1 = max(centerX - floor(width / 2), 1);
    x2 = min(centerX + ceil(width / 2) - 1, imgWidth);
    y1 = max(centerY - floor(height / 2), 1);
    y2 = min(centerY + ceil(height / 2) - 1, imgHeight);

    cropped = image(y1:y2, x1:x2, :);
end
