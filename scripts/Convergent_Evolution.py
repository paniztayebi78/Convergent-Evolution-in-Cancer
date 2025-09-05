import numpy as np
import matplotlib.pyplot as plt
from collections import Counter
from scipy.stats import entropy, pearsonr
import seaborn as sns

# Set random seed for reproducibility
np.random.seed(42)

# Simulation parameters
POPULATION_SIZE = 200
MUTATION_RATE = 0.01
GENERATIONS_SHORT = 50
GENERATIONS_LONG = 500  # Increased for better separation

class TumorEvolutionModel:
    def __init__(self, genome_size, generations, fitness_type='additive'):
        self.genome_size = genome_size
        self.generations = generations
        self.fitness_type = fitness_type
        self.population_size = POPULATION_SIZE
        self.mutation_rate = MUTATION_RATE
        
    def initialize_population(self):
        """Initialize population with all zeros genotype (no standing variation)"""
        return np.zeros((self.population_size, self.genome_size), dtype=int)
    
    def calculate_fitness(self, genotype):
        """Calculate fitness based on genotype"""
        if self.fitness_type == 'additive':
            # Simple additive fitness: more 1s = higher fitness
            return np.sum(genotype)
        elif self.fitness_type == 'peak':
            # Peak fitness at intermediate values (convergent evolution model)
            # Fitness is still sum of bits, but with peak at intermediate value
            ones_count = np.sum(genotype)
            target = self.genome_size // 2  # Peak at 6 1s for 12-bit genome
            # Fitness decreases as you move away from target (6)
            # This allows multiple genotypes with same number of 1s to have same fitness
            distance_from_target = abs(ones_count - target)
            return max(0, target - distance_from_target)  # Peak fitness = 6, decreases with distance
    
    def mutate(self, genotype):
        """Apply random mutations to a genotype"""
        new_genotype = genotype.copy()
        for i in range(len(genotype)):
            if np.random.random() < self.mutation_rate:
                new_genotype[i] = 1 - new_genotype[i]  # flip bit
        return new_genotype
    
    def moran_step(self, population):
        """Single Moran process step: one death, one birth"""
        # Calculate fitness for each individual
        fitnesses = np.array([self.calculate_fitness(ind) for ind in population])
        
        # Avoid division by zero
        if np.sum(fitnesses) == 0:
            birth_probs = np.ones(len(fitnesses)) / len(fitnesses)
        else:
            birth_probs = fitnesses / np.sum(fitnesses)
        
        # Select parent for reproduction (fitness-proportional)
        parent_idx = np.random.choice(len(population), p=birth_probs)
        
        # Select individual to die (random death)
        death_idx = np.random.choice(len(population))
        
        # Create offspring with possible mutations
        offspring = self.mutate(population[parent_idx])
        
        # Replace dead individual with offspring
        new_population = population.copy()
        new_population[death_idx] = offspring
        
        return new_population
    
    def calculate_diversity(self, population):
        """Calculate genotypic and phenotypic diversity using Shannon entropy"""
        # Genotypic diversity: order matters (1,0,0) != (0,1,0)
        genotype_tuples = [tuple(ind) for ind in population]
        genotype_counts = Counter(genotype_tuples)
        genotype_probs = np.array(list(genotype_counts.values())) / self.population_size
        genotypic_entropy = entropy(genotype_probs, base=2)
        
        # Phenotypic diversity: based on fitness scores (phenotype = fitness)
        phenotypes = [self.calculate_fitness(ind) for ind in population]
        phenotype_counts = Counter(phenotypes)
        phenotype_probs = np.array(list(phenotype_counts.values())) / self.population_size
        phenotypic_entropy = entropy(phenotype_probs, base=2)
        
        return genotypic_entropy, phenotypic_entropy
    
    def run_simulation(self):
        """Run full Moran process simulation"""
        population = self.initialize_population()
        
        # Track diversity over time
        genotypic_diversity = []
        phenotypic_diversity = []
        fitness_distributions = []
        
        # Calculate initial diversity
        gd, pd = self.calculate_diversity(population)
        genotypic_diversity.append(gd)
        phenotypic_diversity.append(pd)
        
        # Store fitness distribution for analysis
        fitnesses = [self.calculate_fitness(ind) for ind in population]
        fitness_distributions.append(Counter(fitnesses))
        
        # Run Moran process for specified generations
        # Each generation consists of multiple Moran steps
        steps_per_generation = self.population_size
        
        for gen in range(self.generations):
            # Multiple Moran steps per generation
            for _ in range(steps_per_generation):
                population = self.moran_step(population)
            
            # Calculate diversity after each generation
            gd, pd = self.calculate_diversity(population)
            genotypic_diversity.append(gd)
            phenotypic_diversity.append(pd)
            
            # Store fitness distribution
            fitnesses = [self.calculate_fitness(ind) for ind in population]
            fitness_distributions.append(Counter(fitnesses))
        
        return genotypic_diversity, phenotypic_diversity, fitness_distributions

def plot_oncogene_addiction_models():
    """Plot results for oncogene addiction model (additive fitness)"""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))
    
    # Model 1: 3-bit genome, 50 generations
    model_3bit = TumorEvolutionModel(3, GENERATIONS_SHORT, 'additive')
    gd_3bit, pd_3bit, _ = model_3bit.run_simulation()
    
    axes[0].plot(gd_3bit, label='Genotypic Diversity', color='blue', linewidth=2)
    axes[0].plot(pd_3bit, label='Phenotypic Diversity', color='red', linewidth=2)
    axes[0].set_title('3-bit Genome (50 generations)', fontsize=12)
    axes[0].set_xlabel('Generation')
    axes[0].set_ylabel('Diversity (bits)')
    axes[0].legend()
    axes[0].grid(True, alpha=0.3)
    
    # Model 2: 12-bit genome, 50 generations
    model_12bit_short = TumorEvolutionModel(12, GENERATIONS_SHORT, 'additive')
    gd_12bit_short, pd_12bit_short, _ = model_12bit_short.run_simulation()
    
    axes[1].plot(gd_12bit_short, label='Genotypic Diversity', color='blue', linewidth=2)
    axes[1].plot(pd_12bit_short, label='Phenotypic Diversity', color='red', linewidth=2)
    axes[1].set_title('12-bit Genome (50 generations)', fontsize=12)
    axes[1].set_xlabel('Generation')
    axes[1].set_ylabel('Diversity (bits)')
    axes[1].legend()
    axes[1].grid(True, alpha=0.3)
    
    # Model 3: 12-bit genome, 200 generations
    model_12bit_long = TumorEvolutionModel(12, GENERATIONS_LONG, 'additive')
    gd_12bit_long, pd_12bit_long, fitness_dist = model_12bit_long.run_simulation()
    
    axes[2].plot(gd_12bit_long, label='Genotypic Diversity', color='blue', linewidth=2)
    axes[2].plot(pd_12bit_long, label='Phenotypic Diversity', color='red', linewidth=2)
    axes[2].set_title('12-bit Genome (200 generations)', fontsize=12)
    axes[2].set_xlabel('Generation')
    axes[2].set_ylabel('Diversity (bits)')
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)
    
    plt.suptitle('Oncogene Addiction Model: Genotypic and Phenotypic Diversity', fontsize=14)
    plt.tight_layout()
    plt.show()
    
    return (gd_3bit, pd_3bit), (gd_12bit_short, pd_12bit_short), (gd_12bit_long, pd_12bit_long)

def plot_convergent_evolution_model():
    """Plot results for convergent evolution model (peak fitness)"""
    model_convergent = TumorEvolutionModel(12, GENERATIONS_LONG, 'peak')
    gd_conv, pd_conv, fitness_dist = model_convergent.run_simulation()
    
    plt.figure(figsize=(10, 6))
    
    plt.plot(gd_conv, label='Genotypic Diversity', color='blue', linewidth=2)
    plt.plot(pd_conv, label='Phenotypic Diversity', color='red', linewidth=2)
    plt.title('Convergent Evolution Model\n12-bit Genome (Peak Fitness at 6 ones)', fontsize=12)
    plt.xlabel('Generation')
    plt.ylabel('Diversity (bits)')
    plt.legend()
    plt.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.show()
    
    return gd_conv, pd_conv, fitness_dist

def calculate_statistics(results_dict):
    """Calculate correlation coefficients and other statistics"""
    print("\n=== STATISTICAL SUMMARIES ===")
    
    for model_name, (gd, pd) in results_dict.items():
        correlation, p_value = pearsonr(gd, pd)
        print(f"\n{model_name}:")
        print(f"  Correlation coefficient (Genotypic vs Phenotypic): {correlation:.4f}")
        print(f"  P-value: {p_value:.4f}")
        print(f"  Final genotypic diversity: {gd[-1]:.4f} bits")
        print(f"  Final phenotypic diversity: {pd[-1]:.4f} bits")
        
        # Calculate coupling strength
        if correlation > 0.8:
            coupling = "Strongly Coupled"
        elif correlation > 0.5:
            coupling = "Moderately Coupled"
        elif correlation > 0.2:
            coupling = "Weakly Coupled"
        else:
            coupling = "Decoupled"
        print(f"  Diversity relationship: {coupling}")

def plot_fitness_frequency_over_time(fitness_distributions, title):
    """Plot heatmap of fitness frequencies over time"""
    # Removed - as requested
    pass

# Run all simulations
print("Running Oncogene Addiction Models...")
oncogene_results = plot_oncogene_addiction_models()

print("\nRunning Convergent Evolution Model...")
conv_gd, conv_pd, conv_fitness_dist = plot_convergent_evolution_model()

# Organize results for statistical analysis
results_dict = {
    "3-bit Oncogene Model": oncogene_results[0],
    "12-bit Oncogene Model (50 gen)": oncogene_results[1],
    "12-bit Oncogene Model (200 gen)": oncogene_results[2],
    "12-bit Convergent Evolution Model": (conv_gd, conv_pd)
}

# Calculate and display statistics
calculate_statistics(results_dict)

# Plot fitness frequency heatmaps for convergent evolution model
print("\nGenerating fitness frequency analysis...")
plot_fitness_frequency_over_time(conv_fitness_dist, "Convergent Evolution Model")

# Create comparison plot between oncogene addiction and convergent evolution
plt.figure(figsize=(14, 6))

plt.subplot(1, 2, 1)
plt.plot(oncogene_results[2][0], label='Genotypic', color='blue', linewidth=2)
plt.plot(oncogene_results[2][1], label='Phenotypic', color='red', linewidth=2)
plt.title('Oncogene Addiction Model\n(Additive Fitness)', fontsize=12)
plt.xlabel('Generation')
plt.ylabel('Diversity (bits)')
plt.legend()
plt.grid(True, alpha=0.3)

plt.subplot(1, 2, 2)
plt.plot(conv_gd, label='Genotypic', color='blue', linewidth=2)
plt.plot(conv_pd, label='Phenotypic', color='red', linewidth=2)
plt.title('Convergent Evolution Model\n(Peak Fitness)', fontsize=12)
plt.xlabel('Generation')
plt.ylabel('Diversity (bits)')
plt.legend()
plt.grid(True, alpha=0.3)

plt.suptitle('Comparison: Oncogene Addiction vs Convergent Evolution Models', fontsize=14)
plt.tight_layout()
plt.show()

print("\n=== INTERPRETATION ===")
print("Oncogene Addiction Model: Genotypic and phenotypic diversity should be strongly coupled")
print("Convergent Evolution Model: Genotypic and phenotypic diversity should be decoupled")
print("(High genotypic diversity with low phenotypic diversity indicates convergent evolution)")