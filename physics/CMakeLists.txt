set(BehaviourSource diffusion.cpp dual_behaviour.cpp fracturecriterion.cpp isotropiclineardamage.cpp kelvinvoight.cpp laplacian.cpp lineardamage.cpp linearstiffnessgradient.cpp maxstrain.cpp mohrcoulomb.cpp physics_base.cpp physics.cpp  radialstiffnessgradient.cpp ruptureenergy.cpp stiffness_and_fracture.cpp  stiffness.cpp stiffness_with_imposed_deformation.cpp stiffness_with_variable_imposed_deformation_and_fracture.cpp stiffness_with_variable_imposed_deformation.cpp triple_behaviour.cpp void_form.cpp vonmises.cpp  wave.cpp  weibull_distributed_stiffness.cpp  weibull_distributed_stiffness_with_variable_imposed_deformation_and_fracture.cpp diffusion.h dual_behaviour.h fracturecriterion.h isotropiclineardamage.h kelvinvoight.h laplacian.h lineardamage.h linearstiffnessgradient.h maxstrain.h mohrcoulomb.h physics_base.h physics.h  radialstiffnessgradient.h ruptureenergy.h  stiffness_and_fracture.h  stiffness.h stiffness_with_imposed_deformation.h stiffness_with_variable_imposed_deformation_and_fracture.h stiffness_with_variable_imposed_deformation.h triple_behaviour.h void_form.h vonmises.h  wave.h  weibull_distributed_stiffness.h  weibull_distributed_stiffness_with_variable_imposed_deformation_and_fracture.h)
add_library (AmieBehaviour ${BehaviourSource})







