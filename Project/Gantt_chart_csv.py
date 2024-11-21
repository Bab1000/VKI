import pandas as pd
from datetime import datetime, timedelta

# Define the start and end dates
project_start_date = datetime(2024, 11, 14)
project_end_date = datetime(2025, 6, 26)

# Define the tasks and approximate durations (in days) for each part based on WBS
tasks = [
    # Part 1
    ("Literature Review and Background Analysis", 20),
    ("General Theory of Beta Estimation and Boundary Layer Analysis", 7),
    ("Review of Previous Studies on ICP, BL, Stagline, and OpenFOAM", 7),
    ("Evaluation of Existing Models and Tools", 10),
    ("ICP and BL Code: Applications, Limitations, and Validation", 7),
    ("Stagline Code: Capabilities, Dependencies, and Setup", 7),
    ("OpenFOAM Investigation", 7),
    ("Summary of Insights for Data Comparison", 3),
    ("Tool Comparison and Beta Evaluation", 20),
    ("Computational Beta Generation", 10),
    ("Generate Beta with Stagline Code using ICP I.C.s", 5),
    ("Integrate OpenFOAM for Beta Evaluation (optional)", 5),
    ("Comparative Analysis", 10),
    ("Analysis of ICP, ICP + BL, Stagline, and OpenFOAM Outputs", 7),
    ("Corrections for Velocity Gradient Computation", 3),
    ("Experimental Conditions Expansion", 15),
    ("Experimental Data Handling", 10),
    ("Pre-process Experimental Data", 5),
    ("Reliability Assessment and Beta Validation", 15),
    ("Compare Experimental and Computational Beta Outputs", 7),
    ("Reliability and Consistency Check of Computational Tools", 5),
    ("Identify Reliable Tool(s) for Beta Computation", 3),

    # Part 2
    ("Database Preparation", 20),
    ("Database Structuring", 10),
    ("Define Structure Based on Test Conditions", 5),
    ("Outline Data Requirements", 5),
    ("Data Population", 10),
    ("Populate Database with Simulation Results", 5),
    ("Integrate Experimental Data for Dataset", 5),
    ("Data Consistency Check", 5),
    ("Ensure Consistency Across Data Entries", 5),
    ("Neural Network Model Development", 30),
    ("NN Architecture Definition", 10),
    ("Select NN Structure, Layers, Activation Functions", 5),
    ("Define Input Parameters and Target Outputs", 5),
    ("Training and Optimization", 10),
    ("Train NN Using Database", 5),
    ("Optimize Training Parameters", 5),
    ("Validation and Testing", 10),
    ("Validate NN Against Data", 5),
    ("Assess Generalization and Error Rates", 5),
    
    # Final Phase
    ("Documentation and Report Compilation", 15),
    ("Compile Report", 10),
    ("Presentation Preparation", 10),
    ("Develop Presentation Slides", 5),
    ("Conduct Rehearsal", 5)
]

# Calculate start and end dates for each task and store in DataFrame
task_data = []
current_start_date = project_start_date

for task_name, duration in tasks:
    end_date = current_start_date + timedelta(days=duration)
    task_data.append([task_name, current_start_date, end_date, duration])
    current_start_date = end_date + timedelta(days=1)  # Add a day as buffer

# Create DataFrame
df_gantt = pd.DataFrame(task_data, columns=["Task", "Start Date", "End Date", "Duration (days)"], )

# Save to Excel file
file_path = "/home/jpe/VKI/Project/Gantt_Chart_Project.csv"
df_gantt.to_csv(file_path, index=False, sep=';')

