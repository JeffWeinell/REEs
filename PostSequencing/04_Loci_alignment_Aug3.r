##########################
##### Loads packages #####
##########################

packages.to.load  <- c("ape","seqinr","seqinr","stringr","GenomicRanges","Biostrings","Rsamtools")
packages.to.load2 <- c("rMSA","data.table")
invisible(lapply(packages.to.load, FUN=library, character.only = TRUE))
invisible(lapply(packages.to.load2, FUN=library, character.only = TRUE,lib.loc=packages.dir))
options(stringsAsFactors = FALSE)

###########################
##### Loads functions #####
###########################

source("/Users/Jeff/Google Drive/KU/ExonCapture_LociSelection/SnakeCap_functions.R") ### loads in some necessary functions

####################################
## Step 1 pull out loci and align ##
####################################

###### PARAMETER SETUP ####
threads    <- "10"
# min.taxa <- 4  #min number taxa to keep an alignment # 3 in original version
# min.taxa <- 20 #min number taxa to keep an alignment temporarily changed to 20 (i.e., all individuals must be present)

#Cluster dirs
work.dir         <- "/Volumes/MyPassport/SequenceCapture"
probe.file       <- "/Volumes/MyPassport/SequenceCapture/Weinell_ProbeSet_Snakes_Final_20020Probes_7Oct2018.fasta"
# out.dir        <- "/Volumes/MyPassport/SequenceCapture/Results/Alignments"
out.dir          <- "/Volumes/MyPassport/SequenceCapture/Results/Alignments_longExons_NoUCEs_includeOphiophagus"
species.loci     <- "/Volumes/MyPassport/SequenceCapture/Results/Step3_contigs/SnakeCap1_contigs.fa"        ### this was generated during Step 3
target.loci      <- "/Volumes/MyPassport/SequenceCapture/Weinell_TargetLoci_Snakes_Final_18April2019.txt"   ### JLW added ###
Ophiophagus.loci <- "/Volumes/MyPassport/SequenceCapture/Results/Ophiophagus-hannah_target-loci_30April2019.fa"

UCE.locus.names <- c("WeinellEntry2153-","WeinellEntry2154-","WeinellEntry2155-","WeinellEntry2156-","WeinellEntry2157-","WeinellEntry2158-","WeinellEntry2159-","WeinellEntry2160-","WeinellEntry2161-","WeinellEntry2162-","WeinellEntry2164-","WeinellEntry2165-","WeinellEntry2166-","WeinellEntry2167-","WeinellEntry2168-","WeinellEntry2169-","WeinellEntry2170-","WeinellEntry2171-","WeinellEntry2173-","WeinellEntry2174-","WeinellEntry2175-","WeinellEntry2176-","WeinellEntry2177-","WeinellEntry2178-","WeinellEntry2179-","WeinellEntry2180-","WeinellEntry2181-","WeinellEntry2182-","WeinellEntry2184-","WeinellEntry2185-","WeinellEntry2186-","WeinellEntry2187-","WeinellEntry2188-","WeinellEntry2190-","WeinellEntry2191-","WeinellEntry2194-","WeinellEntry2195-","WeinellEntry2196-","WeinellEntry2197-","WeinellEntry2198-","WeinellEntry2199-","WeinellEntry2201-","WeinellEntry2202-","WeinellEntry2203-","WeinellEntry2204-","WeinellEntry2205-","WeinellEntry2206-","WeinellEntry2207-","WeinellEntry2208-","WeinellEntry2209-","WeinellEntry2210-","WeinellEntry2211-","WeinellEntry2212-","WeinellEntry2213-","WeinellEntry2214-","WeinellEntry2215-","WeinellEntry2217-","WeinellEntry2218-","WeinellEntry2219-","WeinellEntry2220-","WeinellEntry2221-","WeinellEntry2222-","WeinellEntry2223-","WeinellEntry2224-","WeinellEntry2225-","WeinellEntry2226-","WeinellEntry2227-","WeinellEntry2228-","WeinellEntry2229-","WeinellEntry2231-","WeinellEntry2232-","WeinellEntry2233-","WeinellEntry2234-","WeinellEntry2235-","WeinellEntry2236-","WeinellEntry2237-","WeinellEntry2238-","WeinellEntry2239-","WeinellEntry2240-","WeinellEntry2241-","WeinellEntry2242-","WeinellEntry2244-","WeinellEntry2245-","WeinellEntry2246-","WeinellEntry2247-","WeinellEntry2248-","WeinellEntry2249-","WeinellEntry2250-","WeinellEntry2251-","WeinellEntry2252-","WeinellEntry2253-","WeinellEntry2254-","WeinellEntry2255-","WeinellEntry2256-","WeinellEntry2257-","WeinellEntry2258-","WeinellEntry2259-","WeinellEntry2260-","WeinellEntry2261-","WeinellEntry2262-","WeinellEntry2263-","WeinellEntry2264-","WeinellEntry2265-","WeinellEntry2266-","WeinellEntry2267-","WeinellEntry2268-","WeinellEntry2269-","WeinellEntry2270-","WeinellEntry2271-","WeinellEntry2273-","WeinellEntry2274-","WeinellEntry2275-","WeinellEntry2277-","WeinellEntry2278-","WeinellEntry2279-","WeinellEntry2280-","WeinellEntry2281-","WeinellEntry2282-","WeinellEntry2283-","WeinellEntry2284-","WeinellEntry2285-","WeinellEntry2286-","WeinellEntry2287-","WeinellEntry2288-","WeinellEntry2289-","WeinellEntry2290-","WeinellEntry2291-","WeinellEntry2292-","WeinellEntry2293-","WeinellEntry2294-","WeinellEntry2295-","WeinellEntry2297-","WeinellEntry2298-","WeinellEntry2299-","WeinellEntry2300-","WeinellEntry2301-","WeinellEntry2302-","WeinellEntry2303-","WeinellEntry2304-","WeinellEntry2305-","WeinellEntry2306-","WeinellEntry2307-","WeinellEntry2308-","WeinellEntry2309-","WeinellEntry2310-","WeinellEntry2311-","WeinellEntry2312-","WeinellEntry2313-","WeinellEntry2314-","WeinellEntry2315-","WeinellEntry2316-","WeinellEntry2317-","WeinellEntry2318-","WeinellEntry2319-","WeinellEntry2320-","WeinellEntry2321-","WeinellEntry2323-","WeinellEntry2324-","WeinellEntry2325-","WeinellEntry2326-","WeinellEntry2327-","WeinellEntry2328-","WeinellEntry2329-","WeinellEntry2330-","WeinellEntry2331-","WeinellEntry2333-","WeinellEntry2334-","WeinellEntry2335-","WeinellEntry2336-","WeinellEntry2338-","WeinellEntry2339-","WeinellEntry2340-","WeinellEntry2341-","WeinellEntry2342-","WeinellEntry2343-","WeinellEntry2344-","WeinellEntry2345-","WeinellEntry2347-","WeinellEntry2348-","WeinellEntry2349-","WeinellEntry2350-","WeinellEntry2351-","WeinellEntry2352-","WeinellEntry2353-","WeinellEntry2354-","WeinellEntry2355-","WeinellEntry2356-","WeinellEntry2357-","WeinellEntry2358-","WeinellEntry2359-","WeinellEntry2360-","WeinellEntry2362-","WeinellEntry2363-","WeinellEntry2364-","WeinellEntry2365-","WeinellEntry2366-","WeinellEntry2367-","WeinellEntry2368-","WeinellEntry2369-","WeinellEntry2370-","WeinellEntry2371-","WeinellEntry2372-","WeinellEntry2373-","WeinellEntry2375-","WeinellEntry2376-","WeinellEntry2377-","WeinellEntry2378-","WeinellEntry2379-","WeinellEntry2380-","WeinellEntry2381-","WeinellEntry2382-","WeinellEntry2383-","WeinellEntry2384-","WeinellEntry2385-","WeinellEntry2386-","WeinellEntry2387-","WeinellEntry2388-","WeinellEntry2389-","WeinellEntry2390-","WeinellEntry2391-","WeinellEntry2392-","WeinellEntry2393-","WeinellEntry2394-","WeinellEntry2395-","WeinellEntry2396-","WeinellEntry2397-","WeinellEntry2398-","WeinellEntry2399-","WeinellEntry2400-","WeinellEntry2401-","WeinellEntry2402-","WeinellEntry2403-","WeinellEntry2404-","WeinellEntry2405-","WeinellEntry2406-","WeinellEntry2407-","WeinellEntry2409-","WeinellEntry2410-","WeinellEntry2411-","WeinellEntry2412-","WeinellEntry2414-","WeinellEntry2415-","WeinellEntry2416-","WeinellEntry2417-","WeinellEntry2418-","WeinellEntry2420-","WeinellEntry2421-","WeinellEntry2422-","WeinellEntry2423-","WeinellEntry2424-","WeinellEntry2425-","WeinellEntry2426-","WeinellEntry2427-","WeinellEntry2428-","WeinellEntry2429-","WeinellEntry2430-","WeinellEntry2431-","WeinellEntry2432-","WeinellEntry2433-","WeinellEntry2434-","WeinellEntry2435-","WeinellEntry2436-","WeinellEntry2437-","WeinellEntry2438-","WeinellEntry2439-","WeinellEntry2440-","WeinellEntry2441-","WeinellEntry2442-","WeinellEntry2443-","WeinellEntry2444-","WeinellEntry2445-","WeinellEntry2446-","WeinellEntry2447-","WeinellEntry2448-","WeinellEntry2449-","WeinellEntry2450-","WeinellEntry2451-","WeinellEntry2452-","WeinellEntry2453-","WeinellEntry2454-","WeinellEntry2455-","WeinellEntry2456-","WeinellEntry2457-","WeinellEntry2458-","WeinellEntry2459-","WeinellEntry2460-","WeinellEntry2461-","WeinellEntry2462-","WeinellEntry2463-","WeinellEntry2464-","WeinellEntry2466-","WeinellEntry2467-","WeinellEntry2468-","WeinellEntry2469-","WeinellEntry2470-","WeinellEntry2471-","WeinellEntry2472-","WeinellEntry2473-","WeinellEntry2474-","WeinellEntry2476-","WeinellEntry2477-","WeinellEntry2478-","WeinellEntry2479-","WeinellEntry2480-","WeinellEntry2481-","WeinellEntry2482-","WeinellEntry2483-","WeinellEntry2484-","WeinellEntry2486-","WeinellEntry2487-","WeinellEntry2488-","WeinellEntry2489-","WeinellEntry2490-","WeinellEntry2493-","WeinellEntry2494-","WeinellEntry2498-","WeinellEntry2500-","WeinellEntry2501-","WeinellEntry2502-","WeinellEntry2503-","WeinellEntry2504-","WeinellEntry2505-","WeinellEntry2506-","WeinellEntry2507-","WeinellEntry2508-","WeinellEntry2509-","WeinellEntry2510-","WeinellEntry2511-","WeinellEntry2512-","WeinellEntry2513-","WeinellEntry2514-","WeinellEntry2515-","WeinellEntry2516-","WeinellEntry2519-","WeinellEntry2520-","WeinellEntry2521-","WeinellEntry2522-","WeinellEntry2523-","WeinellEntry2524-","WeinellEntry2525-","WeinellEntry2526-","WeinellEntry2527-","WeinellEntry2528-","WeinellEntry2529-","WeinellEntry2530-","WeinellEntry2531-","WeinellEntry2533-","WeinellEntry2534-","WeinellEntry2535-","WeinellEntry2536-","WeinellEntry2537-","WeinellEntry2539-","WeinellEntry2540-","WeinellEntry2541-","WeinellEntry2542-","WeinellEntry2543-","WeinellEntry2544-","WeinellEntry2545-","WeinellEntry2547-","WeinellEntry2548-","WeinellEntry2549-","WeinellEntry2550-","WeinellEntry2551-","WeinellEntry2552-","WeinellEntry2554-","WeinellEntry2555-","WeinellEntry2556-","WeinellEntry2557-","WeinellEntry2558-","WeinellEntry2559-","WeinellEntry2560-","WeinellEntry2561-","WeinellEntry2564-","WeinellEntry2565-","WeinellEntry2566-","WeinellEntry2567-","WeinellEntry2568-","WeinellEntry2569-","WeinellEntry2570-","WeinellEntry2571-","WeinellEntry2572-","WeinellEntry2573-","WeinellEntry2574-","WeinellEntry2575-","WeinellEntry2577-","WeinellEntry2578-","WeinellEntry2579-","WeinellEntry2580-","WeinellEntry2581-","WeinellEntry2582-","WeinellEntry2583-","WeinellEntry2584-","WeinellEntry2585-","WeinellEntry2586-","WeinellEntry2587-","WeinellEntry2588-","WeinellEntry2589-","WeinellEntry2590-","WeinellEntry2591-","WeinellEntry2592-","WeinellEntry2593-","WeinellEntry2594-","WeinellEntry2595-","WeinellEntry2596-","WeinellEntry2597-","WeinellEntry2598-","WeinellEntry2599-","WeinellEntry2600-","WeinellEntry2601-","WeinellEntry2602-","WeinellEntry2603-","WeinellEntry2604-","WeinellEntry2605-","WeinellEntry2606-","WeinellEntry2607-","WeinellEntry2608-","WeinellEntry2609-","WeinellEntry2610-","WeinellEntry2612-","WeinellEntry2613-","WeinellEntry2614-","WeinellEntry2615-","WeinellEntry2616-","WeinellEntry2617-","WeinellEntry2618-","WeinellEntry2619-","WeinellEntry2620-","WeinellEntry2621-","WeinellEntry2622-","WeinellEntry2623-","WeinellEntry2624-","WeinellEntry2625-","WeinellEntry2626-","WeinellEntry2627-","WeinellEntry2628-","WeinellEntry2629-","WeinellEntry2630-","WeinellEntry2631-","WeinellEntry2632-","WeinellEntry2633-","WeinellEntry2634-","WeinellEntry2635-","WeinellEntry2636-","WeinellEntry2637-","WeinellEntry2639-","WeinellEntry2640-","WeinellEntry2641-","WeinellEntry2642-","WeinellEntry2643-","WeinellEntry2644-","WeinellEntry2645-","WeinellEntry2646-","WeinellEntry2647-","WeinellEntry2648-","WeinellEntry2649-","WeinellEntry2650-","WeinellEntry2651-","WeinellEntry2653-","WeinellEntry2654-","WeinellEntry2655-","WeinellEntry2657-","WeinellEntry2658-","WeinellEntry2660-","WeinellEntry2661-","WeinellEntry2662-","WeinellEntry2663-","WeinellEntry2664-","WeinellEntry2665-","WeinellEntry2666-","WeinellEntry2667-","WeinellEntry2669-","WeinellEntry2670-","WeinellEntry2671-","WeinellEntry2672-","WeinellEntry2673-","WeinellEntry2674-","WeinellEntry2675-","WeinellEntry2676-","WeinellEntry2677-","WeinellEntry2678-","WeinellEntry2679-","WeinellEntry2680-","WeinellEntry2681-","WeinellEntry2682-","WeinellEntry2683-","WeinellEntry2684-","WeinellEntry2685-","WeinellEntry2686-","WeinellEntry2687-","WeinellEntry2688-","WeinellEntry2690-","WeinellEntry2691-","WeinellEntry2692-","WeinellEntry2694-","WeinellEntry2695-","WeinellEntry2696-","WeinellEntry2697-","WeinellEntry2698-","WeinellEntry2699-","WeinellEntry2700-","WeinellEntry2701-","WeinellEntry2702-","WeinellEntry2703-","WeinellEntry2704-","WeinellEntry2705-","WeinellEntry2706-","WeinellEntry2708-","WeinellEntry2709-","WeinellEntry2710-","WeinellEntry2711-","WeinellEntry2712-","WeinellEntry2714-","WeinellEntry2715-","WeinellEntry2718-","WeinellEntry2720-","WeinellEntry2721-","WeinellEntry2722-","WeinellEntry2723-","WeinellEntry2724-","WeinellEntry2725-","WeinellEntry2726-","WeinellEntry2727-","WeinellEntry2728-","WeinellEntry2730-","WeinellEntry2731-","WeinellEntry2732-","WeinellEntry2733-","WeinellEntry2734-","WeinellEntry2736-","WeinellEntry2737-","WeinellEntry2738-","WeinellEntry2739-","WeinellEntry2740-","WeinellEntry2741-","WeinellEntry2743-","WeinellEntry2744-","WeinellEntry2745-","WeinellEntry2746-","WeinellEntry2747-","WeinellEntry2748-","WeinellEntry2749-","WeinellEntry2751-","WeinellEntry2753-","WeinellEntry2754-","WeinellEntry2755-","WeinellEntry2757-","WeinellEntry2758-","WeinellEntry2759-","WeinellEntry2760-","WeinellEntry2761-","WeinellEntry2762-","WeinellEntry2763-","WeinellEntry2764-","WeinellEntry2766-","WeinellEntry2767-","WeinellEntry2768-","WeinellEntry2769-","WeinellEntry2770-","WeinellEntry2771-","WeinellEntry2772-","WeinellEntry2773-","WeinellEntry2774-","WeinellEntry2775-","WeinellEntry2776-","WeinellEntry2777-","WeinellEntry2779-","WeinellEntry2780-","WeinellEntry2781-","WeinellEntry2782-","WeinellEntry2783-","WeinellEntry2784-","WeinellEntry2785-","WeinellEntry2786-","WeinellEntry2787-","WeinellEntry2788-","WeinellEntry2789-","WeinellEntry2790-","WeinellEntry2791-","WeinellEntry2792-","WeinellEntry2793-","WeinellEntry2794-","WeinellEntry2795-","WeinellEntry2796-","WeinellEntry2797-","WeinellEntry2798-","WeinellEntry2799-","WeinellEntry2800-","WeinellEntry2801-","WeinellEntry2802-","WeinellEntry2803-","WeinellEntry2804-","WeinellEntry2806-","WeinellEntry2807-","WeinellEntry2808-","WeinellEntry2809-","WeinellEntry2810-","WeinellEntry2811-","WeinellEntry2812-","WeinellEntry2813-","WeinellEntry2814-","WeinellEntry2815-","WeinellEntry2816-","WeinellEntry2817-","WeinellEntry2819-","WeinellEntry2820-","WeinellEntry2821-","WeinellEntry2822-","WeinellEntry2823-","WeinellEntry2824-","WeinellEntry2826-","WeinellEntry2827-","WeinellEntry2828-","WeinellEntry2829-","WeinellEntry2830-","WeinellEntry2831-","WeinellEntry2832-","WeinellEntry2833-","WeinellEntry2834-","WeinellEntry2835-","WeinellEntry2836-","WeinellEntry2837-","WeinellEntry2838-","WeinellEntry2839-","WeinellEntry2840-","WeinellEntry2841-","WeinellEntry2842-","WeinellEntry2844-","WeinellEntry2845-","WeinellEntry2846-","WeinellEntry2847-","WeinellEntry2848-","WeinellEntry2849-","WeinellEntry2850-","WeinellEntry2851-","WeinellEntry2852-","WeinellEntry2853-","WeinellEntry2854-","WeinellEntry2855-","WeinellEntry2856-","WeinellEntry2857-","WeinellEntry2858-","WeinellEntry2859-","WeinellEntry2860-","WeinellEntry2862-","WeinellEntry2863-","WeinellEntry2864-","WeinellEntry2865-","WeinellEntry2866-","WeinellEntry2867-","WeinellEntry2868-","WeinellEntry2869-","WeinellEntry2870-","WeinellEntry2871-","WeinellEntry2873-","WeinellEntry2874-","WeinellEntry2875-","WeinellEntry2876-","WeinellEntry2877-","WeinellEntry2878-","WeinellEntry2879-","WeinellEntry2880-","WeinellEntry2881-","WeinellEntry2882-","WeinellEntry2883-","WeinellEntry2885-","WeinellEntry2886-","WeinellEntry2887-","WeinellEntry2888-","WeinellEntry2889-","WeinellEntry2890-","WeinellEntry2891-","WeinellEntry2892-","WeinellEntry2894-","WeinellEntry2895-","WeinellEntry2896-","WeinellEntry2897-","WeinellEntry2898-","WeinellEntry2899-","WeinellEntry2900-","WeinellEntry2901-","WeinellEntry2903-","WeinellEntry2904-","WeinellEntry2905-","WeinellEntry2906-","WeinellEntry2907-","WeinellEntry2908-","WeinellEntry2909-","WeinellEntry2910-","WeinellEntry2911-","WeinellEntry2912-","WeinellEntry2913-","WeinellEntry2915-","WeinellEntry2916-","WeinellEntry2917-","WeinellEntry2918-","WeinellEntry2919-","WeinellEntry2921-","WeinellEntry2922-","WeinellEntry2923-","WeinellEntry2924-","WeinellEntry2925-","WeinellEntry2927-","WeinellEntry2928-","WeinellEntry2929-","WeinellEntry2930-","WeinellEntry2931-","WeinellEntry2932-","WeinellEntry2933-","WeinellEntry2934-","WeinellEntry2936-","WeinellEntry2937-","WeinellEntry2938-","WeinellEntry2939-","WeinellEntry2940-","WeinellEntry2941-","WeinellEntry2942-","WeinellEntry2943-","WeinellEntry2944-","WeinellEntry2945-","WeinellEntry2946-","WeinellEntry2947-","WeinellEntry2948-","WeinellEntry2949-","WeinellEntry2950-","WeinellEntry2951-","WeinellEntry2952-","WeinellEntry2953-","WeinellEntry2954-","WeinellEntry2955-","WeinellEntry2956-","WeinellEntry2957-","WeinellEntry2958-","WeinellEntry2959-","WeinellEntry2960-","WeinellEntry2961-","WeinellEntry2962-","WeinellEntry2963-","WeinellEntry2964-","WeinellEntry2965-","WeinellEntry2966-","WeinellEntry2967-","WeinellEntry2968-","WeinellEntry2969-","WeinellEntry2971-","WeinellEntry2972-","WeinellEntry2973-","WeinellEntry2974-","WeinellEntry2975-","WeinellEntry2976-","WeinellEntry2977-","WeinellEntry2978-","WeinellEntry2979-","WeinellEntry2980-","WeinellEntry2981-","WeinellEntry2982-","WeinellEntry2983-","WeinellEntry2984-","WeinellEntry2985-","WeinellEntry2986-","WeinellEntry2987-","WeinellEntry2988-","WeinellEntry2989-","WeinellEntry2990-","WeinellEntry2991-","WeinellEntry2992-","WeinellEntry2993-","WeinellEntry2994-","WeinellEntry2995-","WeinellEntry2996-","WeinellEntry2997-","WeinellEntry2998-","WeinellEntry2999-","WeinellEntry3000-","WeinellEntry3001-","WeinellEntry3002-","WeinellEntry3003-","WeinellEntry3005-","WeinellEntry3006-","WeinellEntry3007-","WeinellEntry3008-","WeinellEntry3009-","WeinellEntry3010-","WeinellEntry3011-","WeinellEntry3012-","WeinellEntry3013-","WeinellEntry3014-","WeinellEntry3015-","WeinellEntry3016-","WeinellEntry3017-","WeinellEntry3018-","WeinellEntry3019-","WeinellEntry3020-","WeinellEntry3021-","WeinellEntry3022-","WeinellEntry3023-","WeinellEntry3024-","WeinellEntry3025-","WeinellEntry3026-","WeinellEntry3027-","WeinellEntry3028-","WeinellEntry3030-","WeinellEntry3031-","WeinellEntry3032-","WeinellEntry3034-","WeinellEntry3035-","WeinellEntry3037-","WeinellEntry3038-","WeinellEntry3039-","WeinellEntry3040-","WeinellEntry3042-","WeinellEntry3043-","WeinellEntry3045-","WeinellEntry3046-","WeinellEntry3047-","WeinellEntry3048-","WeinellEntry3049-","WeinellEntry3050-","WeinellEntry3051-","WeinellEntry3052-","WeinellEntry3053-","WeinellEntry3055-","WeinellEntry3056-","WeinellEntry3057-","WeinellEntry3058-","WeinellEntry3059-","WeinellEntry3060-","WeinellEntry3061-","WeinellEntry3062-","WeinellEntry3063-","WeinellEntry3064-","WeinellEntry3065-","WeinellEntry3066-","WeinellEntry3067-","WeinellEntry3068-","WeinellEntry3069-","WeinellEntry3070-","WeinellEntry3071-","WeinellEntry3072-","WeinellEntry3073-","WeinellEntry3075-","WeinellEntry3076-","WeinellEntry3077-","WeinellEntry3078-","WeinellEntry3080-","WeinellEntry3081-","WeinellEntry3082-","WeinellEntry3083-","WeinellEntry3084-","WeinellEntry3085-","WeinellEntry3086-","WeinellEntry3087-","WeinellEntry3088-","WeinellEntry3089-","WeinellEntry3090-","WeinellEntry3091-","WeinellEntry3092-","WeinellEntry3093-","WeinellEntry3094-","WeinellEntry3095-","WeinellEntry3096-","WeinellEntry3098-","WeinellEntry3099-","WeinellEntry3100-","WeinellEntry3101-","WeinellEntry3102-","WeinellEntry3103-","WeinellEntry3105-","WeinellEntry3106-","WeinellEntry3107-","WeinellEntry3108-","WeinellEntry3109-","WeinellEntry3110-","WeinellEntry3111-","WeinellEntry3112-","WeinellEntry3114-","WeinellEntry3115-","WeinellEntry3116-","WeinellEntry3117-","WeinellEntry3118-","WeinellEntry3119-","WeinellEntry3120-","WeinellEntry3121-","WeinellEntry3122-","WeinellEntry3123-","WeinellEntry3124-","WeinellEntry3125-","WeinellEntry3126-","WeinellEntry3127-","WeinellEntry3128-","WeinellEntry3129-","WeinellEntry3130-","WeinellEntry3131-","WeinellEntry3133-","WeinellEntry3134-","WeinellEntry3135-","WeinellEntry3136-","WeinellEntry3137-","WeinellEntry3138-","WeinellEntry3139-","WeinellEntry3140-","WeinellEntry3141-","WeinellEntry3142-","WeinellEntry3143-","WeinellEntry3144-","WeinellEntry3145-","WeinellEntry3150-","WeinellEntry3152")

#START

##############################################################
### Creates output directories if they don't already exist ###
##############################################################

dir.check.create(out.dir)
dir.check.create(paste0(out.dir, "/target-only_untrimmed")) ### no extra flanking region
dir.check.create(paste0(out.dir, "/all-markers_untrimmed")) ### includes bonus flanking region

#########################
### Loads in the data ###
#########################
setwd(work.dir)

# Load sample file from probe matching step
# all.data <- scanFa(FaFile(species.loci))           # loads up fasta file
new.data         <- scanFa(FaFile(species.loci))
Ophiophagus.data <- scanFa(FaFile(Ophiophagus.loci))
all.data         <- merge.contigs <- append(new.data, Ophiophagus.data)
names(all.data)  <- gsub(pattern="_._",replacement="-",x=names(all.data))   ###| JLW added because "|" was a terrible choice for delimiting names

#Gets locus names
bait.loci <- scanFa(FaFile(probe.file))                                     ### loads up fasta file
names(bait.loci) <- gsub(pattern="_._",replacement="-",x=names(bait.loci))  ###| JLW added because | was a terrible choice for delimiting names

reference.loci        <- scanFa(FaFile(target.loci))                       ###| JLW Added
names(reference.loci) <- paste(names(reference.loci),"-reference",sep="")  ###| JLW Added

# if (length(grep("uce-5k-probes", probe.file)) == 1){
#  names(bait.loci)<-gsub("_.*", "", names(bait.loci))
#  bait.loci<-bait.loci[duplicated(names(bait.loci)) != T]
# }

# locus.names <- unique(names(bait.loci))                                          ###| 
locus.names   <- unique(gsub(pattern="-.*",replacement="-",x=names(bait.loci)))    ###| JLW Changed ###
setwd(out.dir)

#Loops through each locus and writes each species to end of file
for (i in 1:length(locus.names)){
 
  #Match probe names to contig names to acquire data
  # match.data<-all.data[grep(pattern = paste(locus.names[i], "_", sep = ""), x = names(all.data))]    ###|
	match.data <-all.data[grep(pattern = locus.names[i], x = names(all.data))]                         ###| JLW Changed ###
	save.name <- gsub("-","",locus.names[i])
  ##############
  # STEP 1: Skip locus if there are too few taxa
  ##############
  if (length(names(match.data)) < min.taxa){                      #### JLW changed from less than or equal to
    print(paste(save.name, " had too few taxa", sep = ""))
    next
  }
  
  ##############
  # STEP 2: Sets up fasta for aligning
  ##############
  
  # names(match.data) <- gsub(pattern = ".*_\\|_", replacement = "", x = names(match.data))  ###|
  names(match.data) <- gsub(pattern = ".*-", replacement = "", x = names(match.data))        ###| JLW Changed ###
  
  # Gets reference locus
  #ref.locus <- bait.loci[grep(pattern = paste(locus.names[i], "$", sep = ""), x = names(bait.loci))]  ###|
  ref.locus <- reference.loci[grep(pattern = locus.names[i], x = names(reference.loci))]               ###| JLW Changed ###
  
  names(ref.locus) <- paste("Reference_Locus")
  final.loci       <- append(match.data, ref.locus)

  ##############################
  #STEP 3: Runs MAFFT to align #
  ##############################
  
  alignment <- mafft(final.loci,param="--localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123 --thread 6")
  alignment <- as(alignment, "DNAStringSet")
  names(alignment) <- names(final.loci)
#  writeXStringSet(alignment, file=paste(save.name,"_align.fa",sep=""))
  
  #Checks for failed mafft run
  if (length(alignment) == 0){
  	next
  }
  
  reversed <- names(alignment)[grep(pattern = "_R_", names(alignment))]
  if (length(reversed[grep(pattern = "Reference_Locus", reversed)]) == 1){
  	alignment<-reverseComplement(alignment)
  }
  names(alignment) <- gsub(pattern = "_R_", replacement = "", x = names(alignment))
  
  #Gets the divergence to make sure not crazy
  diff      <- pairwise.inf.sites(alignment, "Reference_Locus")
  bad.seqs  <- names(diff)[which(diff >= 0.40)]
  rem.align <- alignment[!names(alignment) %in% bad.seqs]
  
  # Moves onto next loop if there are no good sequences
	if (length(rem.align) < as.numeric(min.taxa)){ 
		#Deletes old files
		#system(paste("rm ", out.dir, "/", locus.names[i], "_align.fa ", sep = ""))
		print(paste(save.name, " had too few taxa", sep = ""))
		next
	}
  
  ### re-align if bad seqs removed
  if (length(bad.seqs) != 0){
    #Runs mafft
    alignment <- mafft(rem.align,param="--localpair --maxiterate 1000 --adjustdirection --quiet --op 3 --ep 0.123 --thread 6")
    alignment <- as(alignment, "DNAStringSet")
	names(alignment) <- names(rem.align)
    reversed  <- names(alignment)[grep(pattern = "_R_", names(alignment))]
    if (length(reversed[grep(pattern = "Reference_Locus", reversed)]) == 1){
    	 alignment <- reverseComplement(alignment)
    }
    names(alignment) <- gsub(pattern = "_R_", replacement = "", x = names(alignment))
  } # end bad.seqs if
  
  ##############
  #STEP 4: Save UCE data as the same for both intron and exon + intron datasets since they don't have this
  ##############
  ## this step not really useful for snake data at the moment (although there are UCEs)
  
#  #Detects from the name
#  if (length(grep(pattern = "uce", x=locus.names[i])) == 1){
#    #writes alignment
#    red.align<-alignment[!names(alignment) %in% "Reference_Locus"]
#    new.align<-strsplit(as.character(red.align), "")
#    mat.align<-lapply(new.align, tolower)
#    m.align<-as.matrix(as.DNAbin(mat.align))
#
#    #readies for saving
#    write.phy(m.align, file=paste(out.dir, "/all-markers_untrimmed/", locus.names[i], ".phy", sep = ""), interleave = F)
#    write.phy(m.align, file=paste(out.dir, "/exon-only_untrimmed/", locus.names[i], ".phy", sep = ""), interleave = F)
#    
#    #Deletes old files
#    system(paste("rm ", out.dir, "/", locus.names[i], "_align.fa ", sep = ""))
#    print(paste(locus.names[i], " UCE saved"))
#    next
#  }
  
  ###########################################################################
  # STEP 4: Make target-only alignments (without bonus flanking sequences) ##
  ###########################################################################
  
  #Removes the edge gaps
  ref.aligned <- as.character(alignment['Reference_Locus'])
  not.gaps    <- str_locate_all(ref.aligned, pattern = "[^-]")[[1]][,1]
  ref.start   <- min(not.gaps)
  ref.finish  <- max(not.gaps)
  temp.intron <- subseq(alignment, ref.start, ref.finish)
  
  #Saves prelim exon file
  save.align  <- temp.intron[!names(temp.intron) %in% "Reference_Locus"]
  new.align   <- strsplit(as.character(save.align), "")
  mat.align   <- lapply(new.align, tolower)
  aligned.set <- as.matrix(as.DNAbin(mat.align))
  
  #readies for saving
  write.phy(aligned.set, file=paste(out.dir, "/target-only_untrimmed/", save.name, ".phy", sep = ""), interleave = F)
  
  #########################################################
  # STEP 5: Make alignments including bonus flanking data #
  #########################################################
  
  #writes alignment
  red.align <- alignment[!names(alignment) %in% "Reference_Locus"]
  new.align <- strsplit(as.character(red.align), "")
  mat.align <- lapply(new.align, tolower)
  m.align   <- as.matrix(as.DNAbin(mat.align))
  
  #readies for saving
  write.phy(m.align, file=paste(out.dir, "/all-markers_untrimmed/", save.name, ".phy", sep = ""), interleave = F)
  
  #Deletes old files
  # system(paste("rm ", out.dir, "/", save.name, "_align.fa ", sep = ""))
  print(paste(save.name, " Locus saved"))

}# end big i loop  


#END SCRIPT

