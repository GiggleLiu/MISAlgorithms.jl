function minx(f, vec)
    local xmin = vec[1]
    fmin = f(xmin)
    for j=2:length(vec)
        @inbounds x = vec[j]
        fmin_ = f(x)
        fmin_ < fmin && (fmin = fmin_; xmin=x)
    end
    return xmin
end
