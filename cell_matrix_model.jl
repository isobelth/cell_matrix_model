parameter_dict = Dict("centre_1_position_x" => 1, "centre_1_position_y" => 2,
                      "centre_2_position_x" => 3, "centre_2_position_y" => 4,
                      "current_radius" => 5, "current_orientation" => 6, "current_deformation" => 7 )

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
    i_cell::Array{Float64, 1} = [centre_1_position_x, centre_1_position_y, centre_1_position_x, centre_1_position_y, current_radius, current_orientation, 0.0]
    return i_cell
end
function migrate_i_cell(i_cell::Array{Float64, 1}, maintained_parameters::Array{Float64, 1})
    direction_to_move = 2*pi*rand()
    amount_to_move = maintained_parameters[maintained_parameters_dict["dr_max"]]*rand()
    i_cell[parameter_dict["centre_1_position_x"]] += amount_to_move*cos(direction_to_move)
    i_cell[parameter_dict["centre_1_position_y"]] += amount_to_move*sin(direction_to_move)
    i_cell[parameter_dict["centre_2_position_x"]] = i_cell[parameter_dict["centre_1_position_x"]]
    i_cell[parameter_dict["centre_2_position_y"]] = i_cell[parameter_dict["centre_1_position_y"]]
    return i_cell
end
function grow_i_cell(i_cell::Array{Float64, 1}, maintained_parameters::Array{Float64, 1})
    i_cell[parameter_dict["current_radius"]] += maintained_parameters[maintained_parameters_dict["dR_max"]]*rand()
    if i_cell[parameter_dict["current_radius"]] >  maintained_parameters[maintained_parameters_dict["R_max"]]
        i_cell[parameter_dict["current_radius"]] = maintained_parameters[maintained_parameters_dict["R_max"]]
    end
    return i_cell
end
###################FUNCTIONS FOR M CELLS#############
function create_m_cell(centre_1_position_x::Float64, centre_1_position_y::Float64, centre_2_position_x::Float64, centre_2_position_y::Float64, current_radius::Float64, current_orientation::Float64, current_deformation::Float64)
    m_cell::Array{Float64, 1} = [centre_1_position_x, centre_1_position_y, centre_2_position_x, centre_2_position_y, current_radius, current_orientation, current_deformation]
    return m_cell
end
function migrate_m_cell(m_cell::Array{Float64, 1}, maintained_parameters::Array{Float64, 1})
    direction_to_move = 2*pi*rand()
    amount_to_move = maintained_parameters[maintained_parameters_dict["dr_max"]]*rand()
    m_cell[parameter_dict["centre_1_position_x"]] += amount_to_move*cos(direction_to_move)
    m_cell[parameter_dict["centre_1_position_y"]] += amount_to_move*sin(direction_to_move)
    m_cell[parameter_dict["centre_2_position_x"]] += amount_to_move*cos(direction_to_move)
    m_cell[parameter_dict["centre_2_position_y"]] += amount_to_move*sin(direction_to_move)
    return m_cell
end
function deform_m_cell(m_cell::Array{Float64, 1}, maintained_parameters::Array{Float64, 1})
    amount_to_deform = maintained_parameters[maintained_parameters_dict["dd_max"]]*rand()
    m_cell[parameter_dict["current_deformation"]] += amount_to_deform
    d::Float64=m_cell[parameter_dict["current_deformation"]]
    R::Float64=m_cell[parameter_dict["current_radius"]]
    x::Float64=d/(2*R)
    max_area::Float64 = pi * maintained_parameters[maintained_parameters_dict["R_max"]] * maintained_parameters[maintained_parameters_dict["R_max"]]
    if x<=1
        new_area::Float64 = round(2*R*R * (pi - acos(d/(2*R))+ 0.5*d* sqrt(4*R*R-d*d)), digits=3)
    else
        new_area::Inf16
    end
    while (abs(new_area-max_area))>0.5
        if new_area - max_area >0.1
            R-=0.01
            new_area = 2*R*R * (pi - acos(d/(2*R))+ 0.5*d* sqrt(4*R*R-d*d))
        else
            R+=0.01
            new_area = 2*R*R * (pi - acos(d/(2*R))+ 0.5*d* sqrt(4*R*R-d*d))
        end
    end
    m_cell[parameter_dict["current_radius"]] = round(R, digits=3)
    m_cell[parameter_dict["centre_1_position_x"]] += (amount_to_deform/2)*cos(m_cell[parameter_dict["current_orientation"]])
    m_cell[parameter_dict["centre_1_position_y"]] += (amount_to_deform/2)*sin(m_cell[parameter_dict["current_orientation"]])
    m_cell[parameter_dict["centre_2_position_x"]] += (amount_to_deform/2)*cos(m_cell[parameter_dict["current_orientation"]])
    m_cell[parameter_dict["centre_2_position_y"]] += (amount_to_deform/2)*sin(m_cell[parameter_dict["current_orientation"]])
    return m_cell
end
function rotate_m_cell(m_cell::Array{Float64, 1}, maintained_parameters::Array{Float64, 1})
    amount_to_rotate::Float64 = maintained_parameters[maintained_parameters_dict["dtheta_max"]]*rand()
    m_cell[parameter_dict["current_orientation"]] += amount_to_rotate

    m_cell[parameter_dict["centre_1_position_x"]] += m_cell[parameter_dict["current_deformation"]]*0.5*cos(m_cell[parameter_dict["current_orientation"]])
    m_cell[parameter_dict["centre_1_position_y"]] += m_cell[parameter_dict["current_deformation"]]*0.5*sin(m_cell[parameter_dict["current_orientation"]])
    m_cell[parameter_dict["centre_2_position_x"]] -= m_cell[parameter_dict["current_deformation"]]*0.5*cos(m_cell[parameter_dict["current_orientation"]])
    m_cell[parameter_dict["centre_2_position_y"]] -= m_cell[parameter_dict["current_deformation"]]*0.5*sin(m_cell[parameter_dict["current_orientation"]])
    return m_cell
end
#############FUNCTIONS TO CHANGE CELL TYPE#######
function replace_i_cell_with_m_cell(i_cell::Array{Float64, 1}, maintained_parameters::Array{Float64, 1})
    new_m_cell = create_m_cell(0.0,0.0,0.0,0.0,0.0,0.0,0.0)
    current_orientation = i_cell[parameter_dict["current_orientation"]]#####can update to divide in preferential directions
    R::Float64 = maintained_parameters[maintained_parameters_dict["R_max"]] ##max radius
    max_area = pi * R* R
    d::Float64 = maintained_parameters[maintained_parameters_dict["initial_deformation_proportion"]]*rand() ##initial displacement

    new_area::Float64 = round(2*R*R * (pi - acos(d/(2*R))+ 0.5*d* sqrt(4*R*R-d*d)), digits=3)
    while (abs(new_area-max_area))>0.5
        if new_area - max_area >0.1
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

    return new_m_cell
end
function replace_m_cell_with_two_i_cells(m_cell::Array{Float64, 1}, maintained_parameters::Array{Float64, 1})
    R = m_cell[parameter_dict["current_radius"]]
    current_orientation = m_cell[parameter_dict["current_orientation"]]
    new_i_cell_1_x_position = m_cell[parameter_dict["centre_1_position_x"]]
    new_i_cell_1_y_position = m_cell[parameter_dict["centre_1_position_y"]]
    new_i_cell_2_x_position = m_cell[parameter_dict["centre_2_position_x"]]
    new_i_cell_2_y_position = m_cell[parameter_dict["centre_2_position_y"]]
    new_i_cell_1::Array{Float64, 1} = create_i_cell(new_i_cell_1_x_position,new_i_cell_1_y_position,R,current_orientation)
    new_i_cell_2::Array{Float64, 1} = create_i_cell(new_i_cell_2_x_position,new_i_cell_2_y_position,R,current_orientation)
    created_cells::Array{Float64} = hcat(new_i_cell_1;new_i_cell_2)
    print("new i cells size", size(new_i_cells))
    return new_i_cells
end
############FUNCTION TO SET UP CLUSTER########
function initialise_cluster(i_cell::Array{Float64,1})
    cluster::Array{Float64, 1} = i_cell
    return cluster
end
function make_a_move(cluster::Array{Float64}, maintained_parameters::Array{Float64})
    index_of_cell_to_choose = rand(1: size(cluster,2)) # includes lower and upper number
    chosen_cell = cluster[:,index_of_cell_to_choose]
    random_number = rand(Float64)

    if (chosen_cell[parameter_dict["centre_1_position_x"]] == chosen_cell[parameter_dict["centre_2_position_x"]]) && (
                                    chosen_cell[parameter_dict["centre_1_position_y"]] == chosen_cell[parameter_dict["centre_2_position_y"]])
        #### I CELL PROTOCOL

        if random_number <= 2000/2501
            chosen_cell = grow_i_cell(chosen_cell, maintained_parameters)
            #print("grow\n")
        else
            chosen_cell = migrate_i_cell(chosen_cell, maintained_parameters)
            #print("migrate")
        end
        #print("current=",chosen_cell[parameter_dict["current_radius"]], "max=",maintained_parameters[maintained_parameters_dict["R_max"]], "\n" )
        if chosen_cell[parameter_dict["current_radius"]] == maintained_parameters[maintained_parameters_dict["R_max"]]
            print("replace I with M")
            chosen_cell = replace_i_cell_with_m_cell(chosen_cell, maintained_parameters)
        end

    else
        #### M CELL PROTOCOL
        if random_number <=1500/2506
            chosen_cell=deform_m_cell(chosen_cell, maintained_parameters)
        elseif random_number <= 2200/2506
            chosen_cell=rotate_m_cell(chosen_cell, maintained_parameters)
        else
            chosen_cell=migrate_m_cell(chosen_cell, maintained_parameters)
        end
        if chosen_cell[parameter_dict["current_deformation"]] > chosen_cell[parameter_dict["current_radius"]]

            print("Replace M with 2I")
            two_i_cells = replace_m_cell_with_two_i_cells(chosen_cell, maintained_parameters)
            print("size of new cells=", size(two_i_cells), "\n")

            # chosen_cell = two_i_cells[:,1]
            # other_cell = two_i_cells[:,2]
            #cluster = vcat(cluster; new_cell)
        end
    end
    cluster[:,index_of_cell_to_choose]=chosen_cell
    return cluster
end


function main()
    maintained_parameters = create_maintained_parameters(5.0, 0.05, 0.05, 0.000001, 0.05, 0.005)
    print(maintained_parameters[maintained_parameters_dict["R_max"]])
    cluster_parameters = create_cluster_parameters(5.0,1.0,1.0,2.0,2.0,5.0)
    starting_cell = create_i_cell(1.0,1.0,4.9, 0.0)
    cluster = initialise_cluster(starting_cell)
    print(cluster)
    i = 0
    while i <1000
        cluster = make_a_move(cluster, maintained_parameters)
        #print(cluster, "\n")

        i += 1
    end
    print(cluster)


end

    ###number of rows =  size(cluster, 2)
    # array::Array{Int64} = [1 2 3; 4 5 6]
    # print(array[2,:])
    # new_i_cell = create_i_cell(1.0,1.0,4.5, 0.0)
    # print(new_i_cell)
    # new_i_cell=migrate_i_cell(new_i_cell, maintained_parameters)
    # print(new_i_cell)




main()
