function G_zi = G_zi(zi)
    % Funzione di selettività G(zi) da Ragno et al. (2023)
    G_zi = zeros(size(zi));
    G_zi(zi < 1.35) = 0.002 .* zi(zi < 1.35).^7.5;
    G_zi(zi >= 1.35) = 14.2 .* (1 - 0.894 ./ sqrt(zi(zi >= 1.35))).^4.5;
end