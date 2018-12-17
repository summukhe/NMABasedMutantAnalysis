CXX = g++
CXXFLAGS = -O3 -std=c++11 -fopenmp -g
LDFLAGS = -ldlib -lblas -llapack

targets = generate_config atom_positional_distribution volume_signature electrostatic_signature occlusion_signature sasa_model_train sasa_signature


all:$(targets)

volume_signature: volume_signature.cc grid.h residue_volume.h trajectory_analysis.h config_reader.h
	$(CXX) $(CXXFLAGS)  $< -o $@ $(LDFLAGS)
	
generate_config: generate_config.cc config_reader.h trajectory.h trajectory_analysis.h
	$(CXX) $(CXXFLAGS)  $< -o $@ $(LDFLAGS)

electrostatic_signature: electrostatic_signature.cc grid.h electrostatics.h trajectory_analysis.h
	$(CXX) $(CXXFLAGS)  $< -o $@ $(LDFLAGS)

occlusion_signature: occlusion_signature.cc grid.h residue_volume.h trajectory_analysis.h
	$(CXX) $(CXXFLAGS)  $< -o $@ $(LDFLAGS)

sasa_model_train: sasa_model_train.cc grid.h trajectory_analysis.h sasa_calculation.h
	$(CXX) $(CXXFLAGS)  $< -o $@ $(LDFLAGS)

sasa_signature: sasa_signature.cc grid.h trajectory_analysis.h sasa_calculation.h trajectory.h
	$(CXX) $(CXXFLAGS)  $< -o $@ $(LDFLAGS)
	
atom_positional_distribution: atom_positional_distribution.cc grid.h trajectory_analysis.h
	$(CXX) $(CXXFLAGS)  $< -o $@ $(LDFLAGS)

clean:
	rm -f $(targets)
