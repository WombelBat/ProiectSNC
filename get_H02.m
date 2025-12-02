function [H] = get_H02(zeta,w)
    s=tf('s');
    H= w^2 / (s^2 + 2 * w * zeta * s + w^2);
end