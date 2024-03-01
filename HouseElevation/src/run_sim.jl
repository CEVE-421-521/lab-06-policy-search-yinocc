using Distributions

"""Helper function for trapezoidal rule"""
function trapz(x, y)
    return sum((x[2:end] - x[1:(end - 1)]) .* (y[2:end] + y[1:(end - 1)])) * 0.5
end

"""
Run the model for a given action and SOW

Expected Annual Damages are computed using the trapezoidal rule
"""
function run_sim(a::Action, sow::SOW, p::ModelParams)

    # first, we calculate the cost of elevating the house
    construction_cost = elevation_cost(p.house, a.Δh_ft)

    # we don't need to recalculate the steps of the trapezoidal integral for each year
    storm_surges_ft = range(
        quantile(sow.surge_dist, 0.0005); stop=quantile(sow.surge_dist, 0.9995), length=150
    )

    eads = map(p.years) do year

        # get the sea level for this year
        slr_ft = sow.slr(year)

        # Compute EAD using trapezoidal rule
        pdf_values = pdf.(sow.surge_dist, storm_surges_ft) # probability of each
        depth_ft_gauge = storm_surges_ft .+ slr_ft # flood at gauge
        depth_ft_house = depth_ft_gauge .- (p.house.height_above_gauge_ft + a.Δh_ft) # flood @ house
        damages_frac = p.house.ddf.(depth_ft_house) ./ 100 # damage
        weighted_damages = damages_frac .* pdf_values # weighted damage
        # Trapezoidal integration of weighted damages
        ead = trapz(storm_surges_ft, weighted_damages) * p.house.value_usd
    end

    years_idx = p.years .- minimum(p.years)
    discount_fracs = (1 - sow.discount_rate) .^ years_idx
    ead_npv = sum(eads .* discount_fracs)
    return -(ead_npv + construction_cost)
end

function run_sim_old(a::Action, sow::SOW, p::ModelParams)

    # first, we calculate the cost of elevating the house
    construction_cost = elevation_cost(p.house, a.Δh_ft)

    # next, we calculate expected annual damages for each year
    # map is just a fancy way of writing a for loop across multiple lines
    eads = map(p.years) do year # equivalent to `for year in years`

        # calculate the sea level for this year
        slr_ft = sow.slr(year)

        # compute EAD using Monte Carlo
        storm_surges_ft = rand(sow.surge_dist, 25_000)
        depth_ft_gauge = storm_surges_ft .+ slr_ft
        depth_ft_house = depth_ft_gauge .- (p.house.height_above_gauge_ft + a.Δh_ft)

        # calculate the expected annual damages
        damages_frac = p.house.ddf.(depth_ft_house) ./ 100 # convert to fraction
        mean(damages_frac) * p.house.value_usd # convert to USD
    end

    # finally, we aggregate the costs and benefits to get the net present value
    years_idx = p.years .- minimum(p.years) # 0, 1, 2, 3, .....
    discount_fracs = (1 - sow.discount_rate) .^ years_idx # 1, 1-r, (1-r)^2, (1-r)^3, .....
    ead_npv = sum(eads .* discount_fracs)
    return -(ead_npv + construction_cost)
end
