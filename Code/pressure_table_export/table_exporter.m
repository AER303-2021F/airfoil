latex_table(p_airfoil, 3, 3)
latex_table(p_rake1, 3, 3)
latex_table(p_rake2, 3, 3)

function table = latex_table(dataset, rounding, power)
    table = latex(sym(round(dataset, rounding) * 10^power));
end

