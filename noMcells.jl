ENV["PLOTS_USE_ATOM_PLOTPANE"] = "false"
using Plots
using PyPlot
using PyCall
@pyimport matplotlib.pyplot as pyplt
using VoronoiDelaunay
sp = pyimport("scipy.spatial")

parameter_dict = Dict("centre_1_position_x" => 1, "centre_1_position_y" => 2,
                      "current_radius" => 3)
maintained_parameters_dict = Dict("R_max" => 1, "dR_max" => 2, "dr_max" => 3)
cluster_parameters_dict = Dict("matrix_stiffness" => 1, "attractor_start_x" => 2, "attractor_start_y" => 3,
                                "attractor_end_x" => 4, "attractor_end_y" => 5, "attractor_strength" => 6)
function create_maintained_parameters(R_max::Float64, dR_max::Float64, dr_max::Float64)
  maintained_parameters::Array{Float64, 1} = [R_max, dR_max, dr_max]
  return maintained_parameters
end
function create_cluster_parameters(matrix_stiffness::Float64, attractor_start_x::Float64, attractor_start_y::Float64, attractor_end_x::Float64, attractor_end_y::Float64, attractor_strength::Float64)
    cluster_parameters::Array{Float64, 1} = [matrix_stiffness, attractor_start_x, attractor_start_y, attractor_end_x, attractor_end_y, attractor_strength]
    return cluster_parameters
end
#############FUNCTIONS FOR I CELLS###################
function create_i_cell(centre_1_position_x::Float64, centre_1_position_y::Float64, current_radius::Float64)
    i_cell::Array{Float64, 1} = [centre_1_position_x, centre_1_position_y, current_radius]
    return i_cell
end
function migrate_i_cell(i_cell::Array{Float64, 1}, maintained_parameters::Array{Float64, 1})
    direction_to_move = 2*pi*rand()
    amount_to_move = maintained_parameters[maintained_parameters_dict["dr_max"]]*rand()
    i_cell[parameter_dict["centre_1_position_x"]] += amount_to_move*cos(direction_to_move)
    i_cell[parameter_dict["centre_1_position_y"]] += amount_to_move*sin(direction_to_move)
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


function replace_Rmax_cell_with_two_Rmin_cells(i_cell::Array{Float64, 1}, maintained_parameters::Array{Float64, 1})
    R = maintained_parameters[maintained_parameters_dict["R_max"]]/sqrt(2)
    random_displacement::Float64 = maintained_parameters[maintained_parameters_dict["dr_max"]]*rand()
    direction_to_move = 2*pi*rand()
    new_i_cell_x_position = i_cell[parameter_dict["centre_1_position_x"]]
    new_i_cell_y_position = i_cell[parameter_dict["centre_1_position_y"]]

    ######################can displace them slightly later

    new_i_cell_1::Array{Float64, 1} = create_i_cell((new_i_cell_x_position + random_displacement*cos(direction_to_move)) , (new_i_cell_y_position + random_displacement*sin(direction_to_move)) , R)
    new_i_cell_2::Array{Float64, 1} = create_i_cell(new_i_cell_x_position,new_i_cell_y_position,R)
    new_i_cells::Array{Float64} = hcat(new_i_cell_1 , new_i_cell_2)
    #print("new i cells size", size(new_i_cells))
    return new_i_cells
end
############FUNCTION TO SET UP CLUSTER########
function initialise_cluster(i_cell::Array{Float64,1})
    cluster::Array{Float64, 1} = i_cell
    return cluster
end
######FUNCTION TO FIND NEAREST NEIGHBOURS#####
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
function julia_delaunay_vertices(triang)
    triang_vertices = triang.vertices
    for i::Int64 in (1:1:size(triang_vertices,1))
        for j::Int64 in (1:1:size(triang_vertices,2))
            triang_vertices[i,j] += 1
        end
    end
    return triang_vertices
end
function make_a_move(cluster::Array{Float64}, maintained_parameters::Array{Float64},  number_of_cells::Int64)
    index_of_cell_to_choose = rand(1: size(cluster,2)) # includes lower and upper number
    chosen_cell = cluster[:,index_of_cell_to_choose]
    random_number = rand(Float64)
    neighbour_dict = [0] #this is what will be returned if we have less than three cells
    new_number_of_cells::Int64 = number_of_cells
    #### I CELL PROTOCOL
    if random_number <= 1/2501
        chosen_cell = grow_i_cell(chosen_cell, maintained_parameters)
    else
        chosen_cell = migrate_i_cell(chosen_cell, maintained_parameters)
    end
    cluster[:,index_of_cell_to_choose]=chosen_cell
    if chosen_cell[parameter_dict["current_radius"]] == maintained_parameters[maintained_parameters_dict["R_max"]]
        print("replace I with M")
        two_i_cells= replace_Rmax_cell_with_two_Rmin_cells(chosen_cell, maintained_parameters)
        cluster[:,index_of_cell_to_choose] = two_i_cells[:,1]
        other_cell = two_i_cells[:,2]
        cluster = hcat(cluster, other_cell)
        new_number_of_cells += 1
    end
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
                dx::Float64 = first_cell[parameter_dict["centre_1_position_x"]] - second_cell[parameter_dict["centre_1_position_x"]]
                dy::Float64 = first_cell[parameter_dict["centre_1_position_y"]] - second_cell[parameter_dict["centre_1_position_y"]]
                r12::Float64 = sqrt((dx*dx)+(dy*dy))

                exp_factor::Float64 = exp(-a*(r12-equilibrium_radius))
                potential_for_pair::Float64 = De*(1-(exp_factor*exp_factor))*(1-(exp_factor*exp_factor))
                total_potential_contributed_by_cell += potential_for_pair
            end
        end
        total_potential += total_potential_contributed_by_cell
    end
    return total_potential
end
function main()
    x_coords::Array{Float64} = [0.0]
    y_coords::Array{Float64} = [0.0]
    points::Array = [0]
    dela=sp.Delaunay
    maintained_parameters = create_maintained_parameters(5.0, 0.05, 0.05)
    cluster_parameters = create_cluster_parameters(5.0,1.0,1.0,2.0,2.0,5.0)
    starting_cell = create_i_cell(1.0,1.0,4.9)
    cluster = initialise_cluster(starting_cell)
    neighbour_dict = Dict
    while size(cluster, 2) < 3
        number_of_cells = size(cluster, 2)
        cluster, new_number_of_cells = make_a_move(cluster, maintained_parameters, number_of_cells)
        if number_of_cells == 2 && new_number_of_cells == 3
            x_coords = cluster[parameter_dict["centre_1_position_x"], :]
            y_coords = cluster[parameter_dict["centre_1_position_y"], :]
            points = hcat(x_coords, y_coords)
            triang = dela(points)
            neighbour_dict = make_a_dictionary_of_neighbours(points, triang)
            print("The neighbour dict is", neighbour_dict, "\n")
        end
    end

    energy_before_making_change::Float64 = morse_potential(cluster, neighbour_dict, 1.0, 1.0)
    i::Int64 = 0
    count_accept_new::Int64 = 0
    count_keep_old::Int64 = 0
    #while size(cluster, 2) >= 3 && size(cluster, 2) < 4
    while i < 50000000
        number_of_cells = size(cluster, 2)
        cluster_after_making_change, new_number_of_cells = make_a_move(cluster, maintained_parameters, number_of_cells)
        #print("There were", number_of_cells, "cells and now there are", new_number_of_cells, "cells\n")
        if number_of_cells != new_number_of_cells
            x_coords = cluster_after_making_change[parameter_dict["centre_1_position_x"], :] #only CoM coords
            y_coords = cluster_after_making_change[parameter_dict["centre_1_position_y"], :] #only CoM coords
            points = hcat(x_coords, y_coords)
            triang = dela(points)
            neighbour_dict = make_a_dictionary_of_neighbours(points, triang)
            print("The neighbour dict is", neighbour_dict, "\n")
        end
        energy_after_making_change::Float64 = morse_potential(cluster_after_making_change, neighbour_dict, 1.0, 1.0)
        change_in_energy::Float64 = energy_after_making_change - energy_before_making_change
        x::Float64 = change_in_energy/ (100000)
        random_number = rand()
        exp_delta_E =  exp(-x)
        if i % 5000 == 0
            print("i=", i, "change in energy =", change_in_energy, "random number = ", random_number, "exp-deltaE=", exp_delta_E, "\n")
        end
        if random_number <= exp_delta_E #change in energy < 0 if new state is lower energy
                #############accept new one
            cluster = cluster_after_making_change
            energy_before_making_change = energy_after_making_change
            if change_in_energy >0
                count_accept_new +=1
            end
        else
            count_keep_old += 1
        end
        i += 1
    end
    print("keep new =", count_accept_new, "keep old", count_keep_old, "\n")


    plt.clf()
    ax=pyplt.gca()
    ax.set_xlim((-50, 50))
    ax.set_ylim((-50, 50))
    ax.set_aspect("equal")

    for k::Int64 in (1:1:size(cluster,2))
        cell_to_draw = cluster[:,k]
        ax.add_artist(draw_an_i_cell(cell_to_draw))

    end

    x_coords = cluster[parameter_dict["centre_1_position_x"], :] #only CoM coords
    y_coords = cluster[parameter_dict["centre_1_position_y"], :] #only CoM coords
    points = hcat(x_coords, y_coords)
    triang = dela(points)
    plt.triplot(points[:,1], points[:,2], triang.simplices, color = :blue)
    plt.plot(points[:,1], points[:,2], lw= 0, marker="o", linestyle="")
    fig = gcf()
    display(fig)

    print(transpose(cluster))
    return cluster
end
main()
