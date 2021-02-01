# Standard import statements
import numpy as np
import pandas as pd
import datetime
from dateutil.relativedelta import relativedelta

# Download latest data
def download_data(metrics):
    
    import requests
    
    verbose = True
    today = datetime.date.today() # Today's date
    
    if(verbose): 
        print('Downloading data')
        print('')

    # Form the link
    link = 'data?areaType=region'
    for metric in metrics:
        link = link+'&metric='+metric
    link = link+'&format=csv'
    
    # Full URL and destination file path
    url = 'https://api.coronavirus.data.gov.uk/v2/'+link
    file = 'data/region_'+today.strftime("%Y-%m-%d")+'.csv'
    
    if(verbose): 
        print('URL: %s' % (url))
        print('')
        print('File: %s' % (file))
        print('')

    req = requests.get(url, allow_redirects=True)
    open(file, 'wb').write(req.content) # Write to disk

    if (verbose):
        print('Metrics downloaded:')
        for metric in metrics:
            print(metric)
        print('')
        print('Download complete')
        print('')

# Read in data, no manipulations apart from renaming and deleting columns
def read_data(infile, metrics):
    
    # Parameters
    verbose = False
    
    # Read data into a pandas data frame
    data = pd.read_csv(infile)

    # Print data to screen
    if (verbose):
        print(type(data))
        print(data)
        print('')
        
    # Convert date column to actual date data type
    data.date = pd.to_datetime(data['date'])

    # Remove unnecessary columns
    data.drop('areaType', inplace=True, axis=1)
    data.drop('areaCode', inplace=True, axis=1) 
    
    # Rename columns
    data.rename(columns={'areaName': 'Region'}, inplace=True, errors="raise")
    data.rename(columns=metrics, inplace=True, errors="raise")

    # Print data to screen again
    if (verbose):
        print(type(data))
        print(data)
        print('')

    # Print specific columns to screen
    # TODO: This is probably wrong
    if (verbose):      
        print(data['date'])
        print('')
        for metric in metrics:
            print(data[metrics[metric]])
            print('')
    
    # Return the massaged pandas data frame
    return data

# Sort
def sort_data(df):
  
    df.sort_values(['Region', 'date'], ascending=[True, False], inplace=True)

# Perform calculations on data (assumes organised in date from high to low for each region)
def calculate_data(regions, df, verbose):
    
    # Parameters
    days_roll = 7

    # Utility function for writing out subset of dataframe with comment
    def data_head(df, comment):
        if (verbose):
            print(comment)
            print(df.head(15))
            print()
    
    # Original
    data_head(df, 'Original data')

    # Sort
    sort_data(df)
    data_head(df, 'Sorted data')

    # Calculate rolling cases and deaths (sum over previous week)
    df['Cases_roll_Mead'] = df.apply(lambda x: df.loc[(df.Region == x.Region) & (df.date <= x.date) & (df.date > x.date+relativedelta(days=-days_roll)), 'Cases'].sum(), axis=1)
    df['Deaths_roll_Mead'] = df.apply(lambda x: df.loc[(df.Region == x.Region) & (df.date <= x.date) & (df.date > x.date+relativedelta(days=-days_roll)), 'Deaths'].sum(), axis=1)
    data_head(df, 'Rolling cases and deaths calculated')

    # Calculate doubling times
    df['Cases_roll_past'] = df.apply(lambda x: df.loc[(df.Region == x.Region) & (df.date == x.date+relativedelta(days=-days_roll)), 'Cases_roll'].sum(), axis=1)
    df['Cases_double'] = days_roll*np.log(2.)/np.log(df['Cases_roll']/df['Cases_roll_past'])
    df.drop('Cases_roll_past', inplace=True, axis=1)
    df['Deaths_roll_past'] = df.apply(lambda x: df.loc[(df.Region == x.Region) & (df.date == x.date+relativedelta(days=-days_roll)), 'Deaths_roll'].sum(), axis=1)
    df['Deaths_double'] = days_roll*np.log(2.)/np.log(df['Deaths_roll']/df['Deaths_roll_past'])
    df.drop('Deaths_roll_past', inplace=True, axis=1)
    data_head(df, 'Doubling times calculated')

# Useful information
def useful_info(regions, data):

    # Parameters
    norm_pop = 100000

    # File info
    print('Today\'s date:', datetime.date.today())
    #print('File date:', file_date)
    latest = data['date'].max().strftime("%Y-%m-%d")
    print('Latest date in file:', latest)
    print()

    # Loop over regions
    for region in regions:
        
        # Calculations for per 100,000 population
        pop = regions[region]
        fac = norm_pop/pop
        
        # Region
        print('Region:', region)

        # Isolate regional data
        df = data.loc[data['Region'] == region]
        df.sort_values(['date'], ascending=[False], inplace=True)

        # Date
        print('Date:', df['date'].iloc[0].strftime("%Y-%m-%d"))
        
        # Daily cases
        daily_cases = df['Cases'].iloc[0]
        norm_daily_cases = daily_cases*fac
        print('Daily new cases: %d, or per 100,000 population: %.1f.' % (daily_cases, norm_daily_cases))
        
        # Weekly cases
        weekly_cases = df['Cases_roll'].iloc[0]
        my_weekly_cases = df['Cases_roll_Mead'].iloc[0]    
        if (weekly_cases != my_weekly_cases):
            raise ValueError('My calculation of weekly cases disagrees with official')
        norm_weekly_cases = weekly_cases*fac
        print('Weekly cases: %d, or per 100,000 population: %.1f.' % (weekly_cases, norm_weekly_cases))
        
        # Daily deaths
        daily_deaths = df['Deaths'].iloc[0]
        norm_daily_deaths = daily_deaths*fac
        print('Daily new deaths: %d, or per 100,000 population: %.1f.' % (daily_deaths, norm_daily_deaths))
        
        # Weekly deaths
        weekly_deaths = df['Deaths_roll'].iloc[0]
        my_weekly_deaths = df['Deaths_roll_Mead'].iloc[0]    
        if (weekly_deaths != my_weekly_deaths):
            raise ValueError('My calculation of weekly deaths disagrees with official')
        norm_weekly_deaths = weekly_deaths*fac
        print('Weekly deaths: %d, or per 100,000 population: %.1f.' % (weekly_deaths, norm_weekly_deaths))
        
        # Cases doubling time
        cases_double = df['Cases_double'].iloc[0]
        if (cases_double > 0):
            print('Cases doubling time [days]: %.1f' % (cases_double))
        else:
            print('Cases halving time [days]: %.1f' % (-cases_double))
        print()

# Plot daily data
def plot_bar_data(data, date, start_date, end_date, regions, outfile=None, pop_norm=True, Nmax=None, plot_type='Square'):

    # Imports
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    import matplotlib.ticker as mticker
    import seaborn as sns
    
    # Parameters
    days_in_roll = 7.
    pop_num = 100000
    bar_width = 0.6
    use_seaborn = False

    ### Figure options ###

    # Seaborn
    #if (use_seaborn):
    #    sns.set()
    #else:
    #    sns.reset_orig
    #    matplotlib.rc_file_defaults()
    if (not use_seaborn):
        sns.reset_orig
        matplotlib.rc_file_defaults()

    # Number of plots
    n = len(regions)

    # Size
    if (n == 9):
        if (plot_type == 'Square'):
            figx = 17.; figy = 13.
        elif (plot_type == 'Long'):
            figx = 17.; figy = 50
        else:
            raise ValueError('plot_type must be either Square or Long')
    elif (n == 1):
        figx = 6.; figy = 4.
    else:
        raise ValueError('Only one or nine regions supported')

    # Cases
    plot_cases = True
    case_bar_color = 'cornflowerblue'
    case_line_color = 'b'
    #case_bar_label = 'Positive tests'
    #case_line_label = 'Rolling positive tests'
    case_line_label = 'Positive tests'

    # Hospitalisations 
    plot_hosp = False
    hosp_fac = 10.
    hosp_bar_color = 'green'
    #hosp_line_color = 'green'
    hosp_bar_label = 'Hospital cases [x%d]' % (int(hosp_fac))
    #hosp_line_label = 'Rolling hospital admissions'  

    # Deaths
    plot_deaths = True
    death_fac = 10.
    death_bar_color = 'indianred'
    death_line_color = 'r'
    #death_bar_label = 'Deaths [times %d]' % (int(death_fac))
    #death_line_label = 'Rolling deaths [times %d]' % (int(death_fac))
    death_line_label = 'Deaths [x%d]' % (int(death_fac))

    # Months
    month_color = 'black'
    month_alpha = 0.10

    # Lockdowns
    lockdown_color = 'red'
    lockdown_alpha = 0.25
    lockdown_lab = 'Lockdown'

    # Special dates
    #special_date_color = 'black'
    relax_color = 'green'
    relax_alpha = lockdown_alpha
    relax_lab = 'Relaxation'

    ### ###

    # Relaxations
    plot_relax = False
    Christmas_date = datetime.date(2020, 12, 25)

    # Lockdowns
    plot_lockdowns = True
    Mar_lockdown_start_date = datetime.date(2020, 3, 23)
    Mar_lockdown_end_date = datetime.date(2020, 7, 5)
    Nov_lockdown_start_date = datetime.date(2020, 11, 5)
    Nov_lockdown_end_date = datetime.date(2020, 12, 2)
    Jan_lockdown_start_date = datetime.date(2021, 1, 5)
    Jan_lockdown_end_date = max(data.date)

    # Months
    plot_months = True
    locator_monthstart = mdates.MonthLocator() # Start of every month
    locator_monthmid = mdates.MonthLocator(bymonthday=15) # Middle of every month
    fmt = mdates.DateFormatter('%b') # Specify the format - %b gives us Jan, Feb...
    
    # Plot
    _, ax = plt.subplots(figsize=(figx, figy), sharex=True)
    #plt.tick_params(axis='x', which='minor', bottom=False)
    
    # Loop over regions
    for i, region in enumerate(regions):
        
        if (pop_norm):
            pop_fac = pop_num/regions[region]
        else:
            pop_fac = 1.
        
        if (n == 9):
            if (plot_type == 'Square'):
                plt.subplot(3, 3, i+1)
            elif (plot_type == 'Long'):
                plt.subplot(9, 1, i+1)
            else:
                raise ValueError('Something went wrong with plot_type')
        elif (n == 1):
            plt.subplot(1, 1, 1)
        else:
            raise ValueError('Only supports either one or nine regions')
          
        # Months shading
        if (plot_months):
            for im in [2, 4, 6, 8, 10, 12]:
                plt.axvspan(datetime.date(2020, im, 1), datetime.date(2020, im, 1)+relativedelta(months=+1), 
                            alpha=month_alpha, 
                            color=month_color)
                plt.axvspan(datetime.date(2021, im, 1), datetime.date(2021, im, 1)+relativedelta(months=+1), 
                            alpha=month_alpha, 
                            color=month_color)

        # Lockdowns
        if (plot_lockdowns):
            plt.axvspan(Mar_lockdown_start_date, Mar_lockdown_end_date, 
                        alpha=lockdown_alpha, 
                        color=lockdown_color, 
                        label=lockdown_lab)
            plt.axvspan(Nov_lockdown_start_date, Nov_lockdown_end_date, 
                        alpha=lockdown_alpha, 
                        color=lockdown_color)
            plt.axvspan(Jan_lockdown_start_date, Jan_lockdown_end_date, 
                        alpha=lockdown_alpha, 
                        color=lockdown_color)

        # Important individual dates
        if (plot_relax):
            plt.axvspan(Christmas_date, Christmas_date+relativedelta(days=+1), 
                        color=relax_color, 
                        alpha=relax_alpha, 
                        label=relax_lab)

        # Plot data
        q = "Region == '%s'" % (region) # Query to isolate regions

        # Cases
        if (plot_cases):
            if (use_seaborn):
                sns.barplot(x=data.query(q)['date'], 
                            y=data.query(q)['Cases']*pop_fac,
                            color=case_bar_color)
            else:
                plt.bar(data.query(q)['date'], 
                        data.query(q)['Cases']*pop_fac,
                        width=bar_width,
                        color=case_bar_color)        
            plt.plot(data.query(q)['date'],
                     data.query(q)['Cases_roll']*pop_fac/days_in_roll,
                     color=case_line_color, 
                     label=case_line_label)

        # Hospitalisations
        if (plot_hosp):
            if (use_seaborn):
                sns.barplot(x=data.query(q)['date'], 
                            y=data.query(q)['Hosp']*hosp_fac*pop_fac,                
                            color=hosp_bar_color)
            else:
                plt.bar(data.query(q)['date'], 
                        data.query(q)['Hospital']*hosp_fac*pop_fac, 
                        width=bar_width,       
                        color=hosp_bar_color,
                        label=hosp_bar_label)

        # Deaths
        if (plot_deaths):
            if (use_seaborn):
                sns.barplot(x=data.query(q)['date'], 
                            y=data.query(q)['Deaths']*death_fac*pop_fac,
                            color=death_bar_color)
            else:
                plt.bar(data.query(q)['date'], 
                        data.query(q)['Deaths']*death_fac*pop_fac,
                        width=bar_width,
                        color=death_bar_color)
            plt.plot(data.query(q).date, 
                     death_fac*data.query(q).Deaths_roll*pop_fac/days_in_roll, 
                     color=death_line_color, 
                     label=death_line_label)

        # Ticks and month arragement on x axis
        X = plt.gca().xaxis
        X.set_major_locator(locator_monthstart)
        X.set_major_formatter(fmt)
        X.set_minor_locator(locator_monthmid)
        X.set_major_formatter(mticker.NullFormatter())
        X.set_minor_formatter(fmt)
        ax.tick_params(axis='x', which='minor', bottom=False)
        plt.xlabel('')

        # Axes limits
        plt.xlim(left=start_date, right=end_date)
        if (Nmax != None): 
            plt.ylim(top=Nmax)
        if (pop_norm):
            plt.ylabel('Number per day per 100,000 population')
        else:           
            plt.ylabel('Total number per day')

        # Finalise
        if (n == 1 and region == 'North East'):
            plt.title(region+'\n%s' %(date.strftime("%Y-%m-%d")), x=0.03, y=0.88, loc='Left', bbox=dict(facecolor='w', edgecolor='k'))
        else:
            plt.title(region, x=0.03, y=0.88, loc='Left', bbox=dict(facecolor='w', edgecolor='k'))
        if (n == 9 and i == 0):
            plt.title(date.strftime("%Y-%m-%d"), x=0.97, y=0.88, loc='Right', bbox=dict(facecolor='w', edgecolor='k'))
        if (n == 9 and i == 1) or (n==1 and region == 'North East'): 
            legend = plt.legend(loc='upper right', framealpha=1.)
            legend.get_frame().set_edgecolor('k')

    if(outfile != None):    
        plt.savefig(outfile, dpi=80)
    plt.show(block = False)   

# Plot daily data
def plot_rolling_data(data, date, start_date, end_date, regions, pop_norm=True, plot_type='Cases'):

    import matplotlib.pyplot as plt
    import matplotlib.dates as mdates
    import matplotlib.ticker as mticker
    import seaborn as sns
    sns.set()
    
    # Parameters
    days_in_roll = 7.
    pop_num = 100000

    ### Figure options ###

    if (plot_type == 'Cases_log' or plot_type == 'Deaths_log'):
        log = True
    else:
        log = False

    if (log):
        plot = plt.semilogy
        if (plot_type == 'Cases_log'):
            ymin = 1.
        elif (plot_type == 'Deaths_log'):
            ymin = 1e-2
        else:
            raise ValueError('Plot_type not recognised')
    else:
        plot = plot.plot
        ymin = 0.

    # Size
    figx = 17.; figy = 6.

    # Months
    month_color = 'black'
    month_alpha = 0.10

    # Lockdowns
    lockdown_color = 'red'
    lockdown_alpha = 0.25
    lockdown_lab = 'Lockdown'

    # Special dates
    relax_color = 'green'
    relax_alpha = lockdown_alpha
    relax_lab = 'Relaxation'

    ### ###

    # Relaxations
    plot_relax = False
    Christmas_date = datetime.date(2020, 12, 25)

    # Lockdowns
    plot_lockdowns = True
    Mar_lockdown_start_date = datetime.date(2020, 3, 23)
    Mar_lockdown_end_date = datetime.date(2020, 7, 5)
    Nov_lockdown_start_date = datetime.date(2020, 11, 5)
    Nov_lockdown_end_date = datetime.date(2020, 12, 2)
    Jan_lockdown_start_date = datetime.date(2021, 1, 5)
    Jan_lockdown_end_date = data.date.iloc[0]

    # Months
    plot_months = True
    locator_monthstart = mdates.MonthLocator() # Start of every month
    locator_monthmid = mdates.MonthLocator(bymonthday=15) # Middle of every month
    fmt = mdates.DateFormatter('%b') # Specify the format - %b gives us Jan, Feb...

    # Plot
    _, ax = plt.subplots(figsize=(figx, figy))

    # Months shading
    if (plot_months):
        for im in [2, 4, 6, 8, 10, 12]:
            plt.axvspan(datetime.date(2020, im, 1), datetime.date(2020, im, 1)+relativedelta(months=+1), 
                        alpha=month_alpha, 
                        color=month_color)
            plt.axvspan(datetime.date(2021, im, 1), datetime.date(2021, im, 1)+relativedelta(months=+1), 
                        alpha=month_alpha, 
                        color=month_color)
            
    # Lockdowns
    if (plot_lockdowns):
        plt.axvspan(Mar_lockdown_start_date, Mar_lockdown_end_date, 
                    alpha=lockdown_alpha, 
                    color=lockdown_color, 
                    label=lockdown_lab)
        plt.axvspan(Nov_lockdown_start_date, Nov_lockdown_end_date, 
                    alpha=lockdown_alpha, 
                    color=lockdown_color)
        plt.axvspan(Jan_lockdown_start_date, Jan_lockdown_end_date, 
                    alpha=lockdown_alpha, 
                    color=lockdown_color)
        
    # Important individual dates
    if (plot_relax):
        plt.axvspan(Christmas_date, Christmas_date+relativedelta(days=+1), 
                    color=relax_color, 
                    alpha=relax_alpha, 
                    label=relax_lab)
        
    # Loop over regions
    for i, region in enumerate(regions):
        
        if (pop_norm):
            pop_fac = pop_num/regions[region]
        else:
            pop_fac = 1.

        # Plot data
        q = "Region == '%s'" % (region) # Query to isolate regions      

        # Cases
        if (plot_type == 'Cases' or plot_type == 'Cases_log'):
            plot(data.query(q)['date'],
                     data.query(q)['Cases_roll']*pop_fac/days_in_roll,
                     color='C{}'.format(i), 
                     label=region)
            plt.title('Daily new positive cases: %s' % (date.strftime("%Y-%m-%d")))

        # Deaths
        elif (plot_type == 'Deaths' or plot_type == 'Deaths_log'):
            plot(data.query(q)['date'], 
                     data.query(q)['Deaths_roll']*pop_fac/days_in_roll, 
                     color='C{}'.format(i), 
                     label=region)
            plt.title('Daily new deaths: %s' % (date.strftime("%Y-%m-%d")))

        # Cases doubling
        elif (plot_type == 'Cases_double'):
            plot.plot(data.query(q)['date'],
                     data.query(q)['Cases_double'],
                     color='C{}'.format(i),
                     ls='-',
                     label=region)
            plt.plot(data.query(q)['date'],
                     -data.query(q)['Cases_double'],
                     color='C{}'.format(i),
                     ls='--')
            plt.title('Cases doubling or halving time: %s' % (date.strftime("%Y-%m-%d")))

        # Deaths doubling
        elif (plot_type == 'Deaths_double'):
            plt.plot(data.query(q)['date'],
                     data.query(q)['Deaths_double'],
                     color='C{}'.format(i), 
                     label=region)
            plt.title('Deaths doubling or halving time: %s' % (date.strftime("%Y-%m-%d")))

        else:
            raise ValueError('plot_type specified incorrectly')

    # Ticks and month arragement on x axis
    X = plt.gca().xaxis
    X.set_major_locator(locator_monthstart)
    X.set_major_formatter(fmt)
    X.set_minor_locator(locator_monthmid)
    #ax.xaxis.set_major_formatter(mticker.NullFormatter())
    #ax.xaxis.set_minor_formatter(fmt)
    X.set_major_formatter(mticker.NullFormatter())
    X.set_minor_formatter(fmt)
    ax.tick_params(axis='x', which='minor', bottom=False)

    # Axes limits
    plt.xlim(left=start_date, right=end_date)
    if (plot_type == 'Cases' or plot_type == 'Deaths' or plot_type == 'Cases_log' or plot_type == 'Deaths_log'):
        if (pop_norm):
            plt.ylabel('Number per day per 100,000 population')       
        else:
            plt.ylabel('Total number per day')
        plt.ylim(bottom=ymin)
    elif (plot_type == 'Cases_double' or plot_type == 'Deaths_double'):
        plt.ylabel('Doubling or halving time in days')
        plt.ylim(bottom=0., top=30.)
    else:
        raise ValueError('plot_type specified incorrectly') 

    # Finalise
    plt.legend(loc='upper left')
    plt.show()