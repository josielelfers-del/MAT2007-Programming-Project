import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

##Data Exploration
def data_exploration(file_name):

    print("--- Dataset Information ---")

    #Load the dataset
    try:
        data = pd.read_csv(file_name)

    except FileNotFoundError:
        print("Error: Dataset file was not found.")
        return None

    except Exception:
        print("Error: Dataset could not be read.")
        return None

    ##Uncomment the following print lines to explore the dataset

    # Display the first few rows of the dataset
    #print(data.head())

    # Display column names, note which columns I need (Antibiotic + Resistant Phenotype)
    #print(data.columns)

    # Display information about the dataset, see if there are missing values in the columns of interest
    #print(data.info())

    # Look at the different antibiotics in the dataset
    #print(data["Antibiotic"].unique())

    # Look at the different resistance phenotypes
    #print(data["Resistant Phenotype"].unique())

    print(f"Loading data from: [{file_name}]")
    print(f"Total observations: {len(data)}")
    print(f"Total Variables: {len(data.columns)}")

    return data


##Data Cleaning
def clean_dataset(data):
    print("\n--- Data Cleaning ---")
    
    # standardize capitalization and remove extra spaces
    data["Resistant Phenotype"] = data["Resistant Phenotype"].str.strip().str.capitalize()

    #check for missing values in the "Resistant Phenotype" column
    total_rows = len(data)
    missing_values = data[data["Resistant Phenotype"].isnull()] #find what rows have missing values 

    print(f"Missing resistance phenotypes: {len(missing_values)}")

    # check for missing antibiotic names
    missing_antibiotics = data[data["Antibiotic"].isnull()]
    print(f"Missing antibiotic names: {len(missing_antibiotics)}")

    #Check which data i will use, only rows with nun-null values in the "Resistant Phenotype" column
    clean_data = data[data["Resistant Phenotype"].notnull() & data["Antibiotic"].notnull()] #only keep rows with non-null values in the "Resistant Phenotype" and "Antibiotic" columns
    print(f"Usable Observations: {len(clean_data)}")

    #Check which antibiotics are lefts after cleaning resistance phenotype data
    antibiotics_clean = clean_data["Antibiotic"].unique()

    #Find which antibiotics were removed due to missing resistance phenotype data
    antibiotics_removed = []
    for antibiotic in data["Antibiotic"].unique():
        if antibiotic not in antibiotics_clean:
            antibiotics_removed.append(antibiotic)
    print(f"Number of antibiotics removed: {len(antibiotics_removed)}")
    #print(f"Antibiotics removed due to missing resistance phenotype data: {antibiotics_removed}")
    print(f"Number of antibiotics with usable phenotype data: {len(antibiotics_clean)}")

    return clean_data, antibiotics_clean



##Data Analysis
def data_analysis(clean_data):

    print("\n--- Data Analysis ---")
    # count the number of type of phenotype for each antibiotic
    phenotypes = ["Resistant", "Nonsusceptible", "Susceptible", "Intermediate", "Susceptible-dose dependent"]  

    for phenotype in phenotypes: 
        phenotype_data = clean_data[clean_data["Resistant Phenotype"] == phenotype]
        print(f"Number of {phenotype} phenotypes: {len(phenotype_data)}")

    return phenotypes

def calculate_resistance_rate(antibiotic_data):
    
    #check if there are no observations --> divide by zero error
    if len(antibiotic_data) == 0:
        return 0, 0, 0

    ##Calculate resistance rates for each antibiotic
    resistant_data = antibiotic_data[antibiotic_data["Resistant Phenotype"] == "Resistant"]
    resistance_rate = (len(resistant_data) / len(antibiotic_data)) * 100

    # Calculate standard error of the resistance rate
    proportion_resistant = resistance_rate / 100
    standard_error = ((proportion_resistant * (1 - proportion_resistant)) / len(antibiotic_data)) ** 0.5
    standard_error = standard_error * 100  # convert to percentage
    
    return resistance_rate, len(resistant_data), standard_error ##returns resistant rate, count of resistant observations, and standard error of the resistance rate

def resistance_rates(clean_data, antibiotics_clean):
    print("\n-- Antibiotic Resistance Rates --")
    #calculate the resistance rate for each antibiotic 
    antibiotic_resistance_rates = []
    for antibiotic in antibiotics_clean: 
        antibiotic_data = clean_data[clean_data["Antibiotic"] == antibiotic]
        resistance_rate, resistant_count, standard_error = calculate_resistance_rate(antibiotic_data)
    
        antibiotic_resistance_rates.append((antibiotic, len(antibiotic_data), resistant_count, resistance_rate, standard_error))

    #sort list from highest to lowest resistance rate
    antibiotic_resistance_rates.sort(key=lambda x: x[3], reverse=True) #key=lambda x: x[3] sorts by the 4th element in the list (resistance rate); reverse=True sorts in descending order
    return antibiotic_resistance_rates



##Filtering

def filter_resistance_rates(antibiotic_resistance_rates):

    #exclude antibiotics with less than 100 observations from the final list of resistance rates
    filtered_resistance_rates = []
    for result in antibiotic_resistance_rates:
        if result[1] >= 100:  # result[1] is the total number of observations
            filtered_resistance_rates.append(result)
    print("Minimum sample size for inclusion in final analysis: 100 observations")
    print(f"Number of antibiotics included in final analysis: {len(filtered_resistance_rates)}\n")

    # uncomment the following line to see the filtered resistance rates
    #print(filtered_resistance_rates)



    #Get the top 10 antibiotics with the highest resistance rates
    top_10_antibiotics = filtered_resistance_rates[:10]

    #format into a dataframe
    top_10_table = pd.DataFrame(top_10_antibiotics, columns=["Antibiotic", "Total Observations", "Resistant", "Resistance Rate (%)", "Standard Error (%)"])

    top_10_table["Resistance Rate (%)"] = top_10_table["Resistance Rate (%)"].round(2)
    top_10_table["Standard Error (%)"] = top_10_table["Standard Error (%)"].round(2)

    print("Top 10 Antibiotics by Resistance Rate:\n")
    print(top_10_table.to_string(index=False))

    print("Phenotype distribution of the top 10 antibiotics:\n")

    return filtered_resistance_rates, top_10_table



def calculate_phenotype_distribution(antibiotic_data, phenotypes):

    # Prevent division by zero if the antibiotic has no observations
    if len(antibiotic_data) == 0:
        empty_percentages = []

        # Add 0% for each phenotype
        for phenotype in phenotypes:
            empty_percentages.append(0)

        return empty_percentages

    # calculate the percentage of each phenotype for the antibiotic
    phenotype_percentages = []

    # Go through each possible resistance phenotype
    for phenotype in phenotypes:

        # Filter the antibiotic data to only observations with this phenotype
        phenotype_data = antibiotic_data[
            antibiotic_data["Resistant Phenotype"] == phenotype]

        # Calculate what percentage of the antibiotic observations have this phenotype
        phenotype_percentage = (len(phenotype_data) / len(antibiotic_data)) * 100

        # Store the calculated percentage
        phenotype_percentages.append(phenotype_percentage)

    # Return the percentages for all phenotypes
    return phenotype_percentages



##Calculate phenotype distributions for the Top 10 antibiotics
def phenotype_distributions(top_10_table, clean_data, phenotypes):

    phenotype_distribution = []
    # Go through each antibiotic in the Top 10 table
    for antibiotic in top_10_table["Antibiotic"]:

        # Select all observations belonging to this antibiotic
        antibiotic_data = clean_data[clean_data["Antibiotic"] == antibiotic]

        #calculate the phenotype percentages
        phenotype_percentages = calculate_phenotype_distribution(antibiotic_data, phenotypes)

        # Add the antibiotic name followed by its phenotype percentages to the results list
        phenotype_distribution.append([antibiotic] + phenotype_percentages)



    # Convert the results into a DataFrame
    phenotype_table = pd.DataFrame(phenotype_distribution, columns=["Antibiotic"] + phenotypes)
    for phenotype in phenotypes:
        phenotype_table[phenotype] = phenotype_table[phenotype].round(2)


    #print results
    print("\n-- Phenotype Distribution --")
    print("Phenotype percentages for Top 10 antibiotics:\n")
    print(phenotype_table.to_string(index=False))

    return phenotype_table



def sample_size_analysis(filtered_resistance_rates):

    print("\n-- Sample size vs. resistance rate --")

    ##analysis of the relationship between sample size and resistance rate
    sample_sizes = []
    resistance_rates = []

    for result in filtered_resistance_rates:
        sample_sizes.append(result[1])
        resistance_rates.append(result[3])
    
    # Calculate correlation between sample size and resistance rate --> pearson correlation
    #  make sure there is enough data to calculate correlation
    if len(sample_sizes) >= 2:
        correlation = np.corrcoef(sample_sizes, resistance_rates)[0, 1]
        print(f"Correlation between sample size and resistance rate: {correlation:.3f}")
    else:
        print("Error: Not enough data to calculate correlation.")
        correlation = None

    return sample_sizes, resistance_rates, correlation


##Plotting

def plot_results(top_10_table, sample_sizes, resistance_rates, correlation):

    #Plot the top 10 antibiotics with the highest resistance rates
    plt.figure(figsize=(12, 6))
    plt.barh(top_10_table["Antibiotic"],top_10_table["Resistance Rate (%)"],xerr=top_10_table["Standard Error (%)"],capsize=4) 
    #xerr --> adds error bars to the plot
    plt.gca().invert_yaxis() #Inverts y-axis

    plt.title("Top 10 Antibiotics with the Highest Resistance Rates")
    plt.ylabel("Antibiotic")
    plt.xlabel("Resistance Rate (%)")

    plt.show()


    #plot the relationship between sample size and resistance rate
    if correlation is None: 
        print("Error: Not enough data to plot sample size vs resistance rate.")
        return
    else:
        plt.figure(figsize=(12, 6))
        plt.scatter(sample_sizes, resistance_rates)

        plt.xlabel("Number of Observations")
        plt.ylabel("Resistance Rate (%)")
        plt.title(f"Sample Size vs Resistance Rate (r = {correlation:.3f})")


        plt.show()

    print("\nAnalysis complete.")




def project():
    file_name = "BVBRC_E.coli_Dataset.csv"
    data = data_exploration(file_name)
    if data is None:
        return
    clean_data, antibiotics_clean = clean_dataset(data)
    phenotypes = data_analysis(clean_data)
    antibiotic_resistance_rates = resistance_rates(clean_data, antibiotics_clean)
    filtered_resistance_rates, top_10_table = filter_resistance_rates(antibiotic_resistance_rates)
    phenotype_table = phenotype_distributions(top_10_table, clean_data, phenotypes)
    sample_sizes, resistance_rate_values, correlation = sample_size_analysis(filtered_resistance_rates)
    plot_results(top_10_table, sample_sizes, resistance_rate_values, correlation)

project()
