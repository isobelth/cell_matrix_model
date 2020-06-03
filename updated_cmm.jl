ENV["PLOTS_USE_ATOM_PLOTPANE"] = "false"
using Plots
using PyPlot
using PyCall
@pyimport matplotlib.pyplot as pyplt
using VoronoiDelaunay
sp = pyimport("scipy.spatial")

parameter_dict = Dict("centre_1_position_x" => 1, "centre_1_position_y" => 2,
                      "centre_2_position_x" => 3, "centre_2_position_y" => 4,
                      "current_radius" => 5, "current_orientation" => 6, "current_deformation" => 7, "CoM_position_x" => 8, "CoM_position_y" => 9 )
maintained_parameters_dict = Dict("R_max" => 1, "dR_max" => 2, "dr_max" => 3,
                                  "dtheta_max" => 4, "dd_max" => 5, "initial_deformation_proportion" => 6)
cluster_parameters_dict = Dict("matrix_stiffness" => 1, "attractor_start_x" => 2, "attractor_start_y" => 3,
                                "attractor_end_x" => 4, "attractor_end_y" => 5, "attractor_strength" => 6)
function create_maintained_parameters(R_max::Float64, dR_max::Float64, dr_max::Float64, dtheta_max::Float64, dd_max::Float64, initial_deformation_proportion::Float64)
  maintained_parameters::Array{Float64, 1} = [R_max, dR_max, dr_max, dtheta_max, dd_max, initial_deformation_proportion]
  return maintained_parameters
end
function create_cluster_parameters(matrix_stiffness::Float64, attractor_start_x::Float64, attractor_start_y::Float64, attractor_end_x::Float64, attractor_end_y::Float64, attractor_strength::Float64)
    cluster_parameters::Array{Float64, 1} = [matrix_stiffness, attractor_start_x, attractor_start_y, attractor_end_x, attractor_end_y, attractor_strength]
    return cluster_parameters
end
#############FUNCTIONS FOR I CELLS###################
function create_i_cell(centre_1_position_x::Float64, centre_1_position_y::Float64, current_radius::Float64, current_orientation::Float64)
    i_cell::Array{Float64, 1} = [centre_1_position_x, centre_1_position_y, centre_1_position_x, centre_1_position_y, current_radius, current_orientation, 0.0, centre_1_position_x, centre_1_position_y]
    return i_cell
end
function migrate_i_cell(i_cell::Array{Float64, 1}, maintained_parameters::Array{Float64, 1})
    direction_to_move = 2*pi*rand()
    amount_to_move = maintained_parameters[maintained_parameters_dict["dr_max"]]*rand()
    i_cell[parameter_dict["centre_1_position_x"]] += amount_to_move*cos(direction_to_move)
    i_cell[parameter_dict["centre_1_position_y"]] += amount_to_move*sin(direction_to_move)
    i_cell[parameter_dict["centre_2_position_x"]] = i_cell[parameter_dict["centre_1_position_x"]]
    i_cell[parameter_dict["centre_2_position_y"]] = i_cell[parameter_dict["centre_1_position_y"]]
    i_cell[parameter_dict["CoM_position_x"]] = i_cell[parameter_dict["centre_1_position_x"]]
    i_cell[parameter_dict["CoM_position_y"]] = i_cell[parameter_dict["centre_1_position_y"]]
    return i_cell
end
function grow_i_cell(i_cell::Array{Float64, 1}, maintained_parameters::Array{Float64, 1})
    i_cell[parameter_dict["current_radius"]] += maintained_parameters[maintained_parameters_dict["dR_max"]]*rand()
    if i_cell[parameter_dict["current_radius"]] >  maintained_parameters[maintained_parameters_dict["R_max"]]
        i_cell[parameter_dict["current_radius"]] = maintained_parameters[maintained_parameters_dict["R_max"]]
    end
    return i_cell
end
function draw_an_i_cell(i_cell::Array{Float64})
    x_pos = i_cell[parameter_dict["centre_1_position_x"]]
    y_pos = i_cell[parameter_dict["centre_1_position_y"]]
    radius= i_cell[parameter_dict["current_radius"]]
    circle = pyplt.Circle((x_pos, y_pos), radius, alpha=0.1 , color=:green, lw=0.3)
    return circle
end
###################FUNCTIONS FOR M CELLS#############
function create_m_cell(centre_1_position_x::Float64, centre_1_position_y::Float64, centre_2_position_x::Float64, centre_2_position_y::Float64, current_radius::Float64, current_orientation::Float64, current_deformation::Float64)
    m_cell::Array{Float64, 1} = [centre_1_position_x, centre_1_position_y, centre_2_position_x, centre_2_position_y, current_radius, current_orientation, current_deformation, 0.5*(centre_1_position_x+centre_2_position_x),0.5*(centre_1_position_y+centre_2_position_y)]
    return m_cell
end
function migrate_m_cell(m_cell::Array{Float64, 1}, maintained_parameters::Array{Float64, 1})
    direction_to_move = 2*pi*rand()
    amount_to_move = maintained_parameters[maintained_parameters_dict["dr_max"]]*rand()
    m_cell[parameter_dict["centre_1_position_x"]] += amount_to_move*cos(direction_to_move)
    m_cell[parameter_dict["centre_1_position_y"]] += amount_to_move*sin(direction_to_move)
    m_cell[parameter_dict["centre_2_position_x"]] += amount_to_move*cos(direction_to_move)
    m_cell[parameter_dict["centre_2_position_y"]] += amount_to_move*sin(direction_to_move)
    m_cell[parameter_dict["CoM_position_x"]]=0.5*(m_cell[parameter_dict["centre_1_position_x"]]+m_cell[parameter_dict["centre_2_position_x"]])
    m_cell[parameter_dict["CoM_position_y"]]=0.5*(m_cell[parameter_dict["centre_1_position_y"]]+m_cell[parameter_dict["centre_2_position_y"]])
    return m_cell
end
function deform_m_cell(m_cell::Array{Float64, 1}, maintained_parameters::Array{Float64, 1})
    amount_to_deform = maintained_parameters[maintained_parameters_dict["dd_max"]]*rand()
    d::Float64 =  m_cell[parameter_dict["current_deformation"]] + amount_to_deform
    print("Deforming by", amount_to_deform, "new deformation =",  m_cell[parameter_dict["current_deformation"]], "\n")
    R::Float64= m_cell[parameter_dict["current_radius"]]
    #x::Float64=d/(2*R)
    max_area::Float64 = round(pi * maintained_parameters[maintained_parameters_dict["R_max"]] * maintained_parameters[maintained_parameters_dict["R_max"]], digits=4)
    if (d/(2*R))<=1
        new_area::Float64 = round(2*R*R * (pi - acos(d/(2*R))+ (0.5*d* sqrt(4*R*R-d*d))), digits=5)
    else
        new_area::Inf16
    end
    while (abs(new_area-max_area))>0.5
        if new_area - max_area >0.5
            R-=0.01
            new_area = round(2*R*R * (pi - acos(d/(2*R))+ (0.5*d* sqrt(4*R*R-d*d))), digits=5)
        else
            R+=0.01
            new_area = round(2*R*R * (pi - acos(d/(2*R))+ (0.5*d* sqrt(4*R*R-d*d))), digits=5)
        end
    end
    m_cell[parameter_dict["current_deformation"]] = d
    m_cell[parameter_dict["current_radius"]] = round(R, digits=4)
    m_cell[parameter_dict["centre_1_position_x"]] += (amount_to_deform/2)*cos(m_cell[parameter_dict["current_orientation"]])
    m_cell[parameter_dict["centre_1_position_y"]] += (amount_to_deform/2)*sin(m_cell[parameter_dict["current_orientation"]])
    m_cell[parameter_dict["centre_2_position_x"]] += (amount_to_deform/2)*cos(m_cell[parameter_dict["current_orientation"]])
    m_cell[parameter_dict["centre_2_position_y"]] += (amount_to_deform/2)*sin(m_cell[parameter_dict["current_orientation"]])
    m_cell[parameter_dict["CoM_position_x"]]=0.5*(m_cell[parameter_dict["centre_1_position_x"]]+m_cell[parameter_dict["centre_2_position_x"]])
    m_cell[parameter_dict["CoM_position_y"]]=0.5*(m_cell[parameter_dict["centre_1_position_y"]]+m_cell[parameter_dict["centre_2_position_y"]])
    return m_cell
end
function rotate_m_cell(m_cell::Array{Float64, 1}, maintained_parameters::Array{Float64, 1})
    amount_to_rotate::Float64 = maintained_parameters[maintained_parameters_dict["dtheta_max"]]*rand()
    m_cell[parameter_dict["current_orientation"]] += amount_to_rotate

    m_cell[parameter_dict["centre_1_position_x"]] += m_cell[parameter_dict["current_deformation"]]*0.5*cos(m_cell[parameter_dict["current_orientation"]])
    m_cell[parameter_dict["centre_1_position_y"]] += m_cell[parameter_dict["current_deformation"]]*0.5*sin(m_cell[parameter_dict["current_orientation"]])
    m_cell[parameter_dict["centre_2_position_x"]] -= m_cell[parameter_dict["current_deformation"]]*0.5*cos(m_cell[parameter_dict["current_orientation"]])
    m_cell[parameter_dict["centre_2_position_y"]] -= m_cell[parameter_dict["current_deformation"]]*0.5*sin(m_cell[parameter_dict["current_orientation"]])
    m_cell[parameter_dict["CoM_position_x"]]=0.5*(m_cell[parameter_dict["centre_1_position_x"]]+m_cell[parameter_dict["centre_2_position_x"]])
    m_cell[parameter_dict["CoM_position_y"]]=0.5*(m_cell[parameter_dict["centre_1_position_y"]]+m_cell[parameter_dict["centre_2_position_y"]])
    return m_cell
end
function draw_an_m_cell(m_cell::Array{Float64})
    x_pos1 = m_cell[parameter_dict["centre_1_position_x"]]
    y_pos1 = m_cell[parameter_dict["centre_1_position_y"]]
    x_pos2 = m_cell[parameter_dict["centre_2_position_x"]]
    y_pos2 = m_cell[parameter_dict["centre_2_position_y"]]
    radius= m_cell[parameter_dict["current_radius"]]
    print("The radius of the m cells to draw is", radius, "\n")
    circle1 = pyplt.Circle((x_pos1, y_pos1), radius, alpha = 0.1 , color = :blue, lw=0.3)
    circle2 = pyplt.Circle((x_pos2, y_pos2), radius, alpha=0.1 , color = :blue, lw=0.3)
    circles::Array{} = [circle1, circle2]
    return circles
end
function replace_i_cell_with_m_cell(i_cell::Array{Float64, 1}, maintained_parameters::Array{Float64, 1})
    new_m_cell = create_m_cell(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    current_orientation = i_cell[parameter_dict["current_orientation"]]#####can update to divide in preferential directions
    R::Float64 = maintained_parameters[maintained_parameters_dict["R_max"]] ##max radius
    max_area = pi * R* R
    d::Float64 = maintained_parameters[maintained_parameters_dict["initial_deformation_proportion"]]*rand() ##initial displacement

    new_area::Float64 = round(2*R*R * (pi - acos(d/(2*R))+ 0.5*d* sqrt(4*R*R-d*d)), digits=3)
    while (abs(new_area-max_area))>0.5
        if new_area - max_area >0.5
            R-=0.01
            new_area = 2*R*R * (pi - acos(d/(2*R))+ 0.5*d* sqrt(4*R*R-d*d))
        else
            R+=0.01
            new_area = 2*R*R * (pi - acos(d/(2*R))+ 0.5*d* sqrt(4*R*R-d*d))
        end
    end
    new_m_cell[parameter_dict["current_radius"]] = round(R, digits=3)
    new_m_cell[parameter_dict["current_orientation"]] = current_orientation
    new_m_cell[parameter_dict["current_deformation"]] = d
    new_m_cell[parameter_dict["centre_1_position_x"]] = i_cell[parameter_dict["centre_1_position_x"]] + (d)*0.5*cos(current_orientation)
    new_m_cell[parameter_dict["centre_1_position_y"]] = i_cell[parameter_dict["centre_1_position_y"]] + (d)*0.5*sin(current_orientation)
    new_m_cell[parameter_dict["centre_2_position_x"]] = i_cell[parameter_dict["centre_1_position_x"]] - (d)*0.5*cos(current_orientation)
    new_m_cell[parameter_dict["centre_2_position_y"]] = i_cell[parameter_dict["centre_1_position_y"]] - (d)*0.5*sin(current_orientation)
    new_m_cell[parameter_dict["CoM_position_x"]]=(0.5*new_m_cell[parameter_dict["centre_1_position_x"]]+new_m_cell[parameter_dict["centre_2_position_x"]])
    new_m_cell[parameter_dict["CoM_position_y"]]=(0.5*new_m_cell[parameter_dict["centre_1_position_y"]]+new_m_cell[parameter_dict["centre_2_position_y"]])

    return new_m_cell
end
function replace_m_cell_with_two_i_cells(m_cell::Array{Float64, 1}, maintained_parameters::Array{Float64, 1})
    r = m_cell[parameter_dict["current_radius"]]
    current_orientation = m_cell[parameter_dict["current_orientation"]]
    new_i_cell_1_x_position = m_cell[parameter_dict["centre_1_position_x"]]
    new_i_cell_1_y_position = m_cell[parameter_dict["centre_1_position_y"]]
    new_i_cell_2_x_position = m_cell[parameter_dict["centre_2_position_x"]]
    new_i_cell_2_y_position = m_cell[parameter_dict["centre_2_position_y"]]
    new_i_cell_1::Array{Float64, 1} = create_i_cell(new_i_cell_1_x_position,new_i_cell_1_y_position,r,current_orientation)
    new_i_cell_2::Array{Float64, 1} = create_i_cell(new_i_cell_2_x_position,new_i_cell_2_y_position,r,current_orientation)
    new_i_cells::Array{Float64} = hcat(new_i_cell_1 , new_i_cell_2)
    return new_i_cells
end
function initialise_cluster(i_cell::Array{Float64,1})
    cluster::Array{Float64, 1} = i_cell
    return cluster
end
function find_neighbours(pindex, triang)
        pindex -=1
        neighbours::Array{Int64} = []
        for i::Int64 in (1:1:size(triang.vertices,1))
                for j::Int64 in (1:1:3)
                        for k::Int64 in (1:1:3)
                                if triang.vertices[i,j] == pindex
                        #if triang.vertices[i,j] != pindex
                                        if j != k
                                                neighbours=vcat(neighbours, triang.vertices[i,k] + 1)
                                        end
                                end
                        end
                end
        end
        neighbours = unique(neighbours)
        return neighbours
end
function make_a_dictionary_of_neighbours(points::Array{Float64}, triang)
    neighbours_dictionary = Dict()
    for i::Int64 in (1:1:size(points,1))
                #key_for_dict = string(i)
            key_for_dict = i
            element_for_dict = find_neighbours(i,triang)
            neighbours_dictionary[key_for_dict] = element_for_dict
    end
    return neighbours_dictionary
end
function make_a_move(cluster::Array{Float64}, maintained_parameters::Array{Float64}, number_of_cells::Int64)
    index_of_cell_to_choose = rand(1: size(cluster,2)) # includes lower and upper number
    chosen_cell = cluster[:,index_of_cell_to_choose]
    random_number = rand(Float64)
    new_number_of_cells::Int64 = number_of_cells
    if (chosen_cell[parameter_dict["centre_1_position_x"]] == chosen_cell[parameter_dict["centre_2_position_x"]]) && (
                                    chosen_cell[parameter_dict["centre_1_position_y"]] == chosen_cell[parameter_dict["centre_2_position_y"]])
        #### I CELL PROTOCOL
        if random_number <= 1/2501
            chosen_cell = grow_i_cell(chosen_cell, maintained_parameters)
        else
            chosen_cell = migrate_i_cell(chosen_cell, maintained_parameters)
        end
        if chosen_cell[parameter_dict["current_radius"]] == maintained_parameters[maintained_parameters_dict["R_max"]]
            print("replace I with M")
            chosen_cell = replace_i_cell_with_m_cell(chosen_cell, maintained_parameters)
        end
        #cluster[:,index_of_cell_to_choose]=chosen_cell
    else
        #### M CELL PROTOCOL
        if random_number <=1/2506
            print("cell index is", index_of_cell_to_choose, "\n")
            chosen_cell=deform_m_cell(chosen_cell, maintained_parameters)

        elseif random_number <= 5/2506
            chosen_cell=rotate_m_cell(chosen_cell, maintained_parameters)
        else
            chosen_cell=migrate_m_cell(chosen_cell, maintained_parameters)
        end
        Rmrt::Float64 = maintained_parameters[maintained_parameters_dict["R_max"]]/sqrt(2)
        #print("current radius = ", chosen_cell[parameter_dict["current_radius"]], "Rmrt=", Rmrt, "\n")
        if chosen_cell[parameter_dict["current_radius"]] <=Rmrt

            print("Replace M with 2I \n")
            two_i_cells = replace_m_cell_with_two_i_cells(chosen_cell, maintained_parameters)
            chosen_cell = two_i_cells[:,1]
            other_cell = two_i_cells[:,2]
            #cluster[:,index_of_cell_to_choose]=chosen_cell
            cluster = hcat(cluster, other_cell)
            new_number_of_cells = number_of_cells + 1
        end
        #cluster[:,index_of_cell_to_choose]=chosen_cell
    end
    cluster[:,index_of_cell_to_choose]=chosen_cell
    return cluster, new_number_of_cells
end
function morse_potential(cluster::Array{Float64}, neighbour_dict::Dict, De::Float64, a::Float64)
    total_potential::Float64 = 0
    for i::Int64 in (1:1:size(cluster,2))
        first_cell = cluster[:, i]
        neighbours_list = neighbour_dict[i]
        total_potential_contributed_by_cell = 0
        for j::Int64 in (1:1:size(neighbours_list,1))
            if i < j
                second_cell = cluster[:, j]
                equilibrium_radius::Float64 = first_cell[parameter_dict["current_radius"]] + second_cell[parameter_dict["current_radius"]]
                dx::Float64 = first_cell[parameter_dict["CoM_position_x"]] - second_cell[parameter_dict["CoM_position_x"]]
                dy::Float64 = first_cell[parameter_dict["CoM_position_y"]] - second_cell[parameter_dict["CoM_position_y"]]
                r12::Float64 = sqrt((dx*dx)+(dy*dy))

                exp_factor::Float64 = exp(-a*(r12-equilibrium_radius))
                potential_for_pair::Float64 = De*(1-(exp_factor*exp_factor))*(1-(exp_factor*exp_factor))
                #print("i=", i, "j=", j, "potential for pair=", potential_for_pair, "\n")
                total_potential_contributed_by_cell += potential_for_pair
                #print("running total for cell ", i, "=", total_potential_contributed_by_cell, "\n")
            end
        end
        #print("i=", i, "total potential from cell =", total_potential_contributed_by_cell, "\n")
        total_potential += total_potential_contributed_by_cell
    end
    #print("total potential from cluster =", total_potential, "\n")
    return total_potential
end
function main()
    x_coords::Array{Float64} = [0.0]
    y_coords::Array{Float64} = [0.0]
    points::Array = [0]
    dela=sp.Delaunay
    maintained_parameters = create_maintained_parameters(5.0, 0.05, 0.05, 0.000001, 0.05, 0.005)
    cluster_parameters = create_cluster_parameters(5.0,1.0,1.0,2.0,2.0,5.0)
    starting_cell = create_i_cell(1.0,1.0,4.9, 0.0)
    cluster = initialise_cluster(starting_cell)
    neighbour_dict = Dict
    while size(cluster, 2) < 3
        number_of_cells = size(cluster, 2)
        cluster, new_number_of_cells = make_a_move(cluster, maintained_parameters, number_of_cells)
        if number_of_cells == 2 && new_number_of_cells == 3
            x_coords = cluster[parameter_dict["CoM_position_x"], :] #only CoM coords
            y_coords = cluster[parameter_dict["CoM_position_y"], :] #only CoM coords
            points = hcat(x_coords, y_coords)
            triang = dela(points)
            neighbour_dict = make_a_dictionary_of_neighbours(points, triang)
            print("The neighbour dict is", neighbour_dict, "\n")
        end
    end

    energy_before_making_change::Float64 = morse_potential(cluster, neighbour_dict, 5.0, 1.0)
    i::Int64 = 0
    count_accept_new::Int64 = 0
    count_keep_old::Int64 = 0
    #while size(cluster, 2) >= 3 && size(cluster, 2) < 4
    while i < 1500
        number_of_cells = size(cluster, 2)
        cluster_after_making_change, new_number_of_cells = make_a_move(cluster, maintained_parameters, number_of_cells)
        print("There were", number_of_cells, "cells and now there are", new_number_of_cells, "cells\n")
        if number_of_cells != new_number_of_cells
            x_coords = cluster_after_making_change[parameter_dict["CoM_position_x"], :] #only CoM coords
            y_coords = cluster_after_making_change[parameter_dict["CoM_position_y"], :] #only CoM coords
            points = hcat(x_coords, y_coords)
            #dela = scipy.spatial.Delaunay
            triang = dela(points)
            neighbour_dict = make_a_dictionary_of_neighbours(points, triang)
            print("The neighbour dict is", neighbour_dict, "\n")
        end
        energy_after_making_change::Float64 = morse_potential(cluster_after_making_change, neighbour_dict, 1.0, 1.0)
        change_in_energy::Float64 = energy_after_making_change - energy_before_making_change
        x::Float64 = change_in_energy/ (0.0001)
        random_number = rand()
        exp_delta_E =  (1 - (x) + (0.5*x*x))
        print("i=", i, "change in energy =", change_in_energy, "random number = ", random_number, "exp-deltaE=", exp_delta_E, "\n")

        if random_number <= exp_delta_E #change in energy < 0 if new state is lower energy
                #############accept new one
            cluster = cluster_after_making_change
            energy_before_making_change = energy_after_making_change
            count_accept_new +=1
        else
            count_keep_old += 1
        end
        i += 1
    end
    print("keep new =", count_accept_new, "keep old", count_keep_old, "\n")


    plt.clf()
    ax=pyplt.gca()
    ax.set_xlim((-30, 30))
    ax.set_ylim((-30, 30))
    ax.set_aspect("equal")

    for k::Int64 in (1:1:size(cluster,2))
        cell_to_draw = cluster[:,k]
        if (cell_to_draw[parameter_dict["centre_1_position_x"]] == cell_to_draw[parameter_dict["centre_2_position_x"]]) && (
                                        cell_to_draw[parameter_dict["centre_1_position_y"]] == cell_to_draw[parameter_dict["centre_2_position_y"]])
            #I CELL PROTOCOL
            ax.add_artist(draw_an_i_cell(cell_to_draw))
        else
            circles = draw_an_m_cell(cell_to_draw)
            current_cell = circles[1]
            ax.add_artist(current_cell)
            current_cell = circles[2]
            ax.add_artist(current_cell)
        end
    end
    triang = dela(points)
    plt.triplot(points[:,1], points[:,2], triang.simplices, color = :blue)
    plt.plot(points[:,1], points[:,2], lw= 0, marker="o", linestyle="")
    fig = gcf()
    display(fig)

    print(transpose(cluster))
end

main()
