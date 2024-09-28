import streamlit as st

# Define a dictionary with detailed disease information
disease_info = {
    "Alzheimer's Disease": """
Alzheimer's disease is a progressive neurological disorder that causes brain cells to degenerate and die. It is the most common cause of dementia, characterized by a continuous decline in thinking, behavioral, and social skills. 

**Symptoms:**
- Memory loss that disrupts daily life
- Difficulty planning and solving problems
- Confusion with time or place
- Difficulty completing familiar tasks at home, at work, or at leisure
- Trouble understanding visual images and spatial relationships
- New problems with words in speaking or writing
- Misplacing things and losing the ability to retrace steps
- Decreased or poor judgment
- Withdrawal from work or social activities
- Changes in mood and personality

**Pathophysiology:**
- Accumulation of amyloid-beta plaques and tau tangles in the brain
- Neuroinflammation and loss of synaptic connections
- Progressive neuronal loss in the hippocampus and other brain regions

**Current Treatments:**
- Cholinesterase inhibitors (e.g., donepezil, rivastigmine) to increase acetylcholine levels
- NMDA receptor antagonist (memantine) to regulate glutamate activity
- Medications to manage symptoms like depression, agitation, and sleep disturbances
- Lifestyle interventions including cognitive stimulation, physical activity, and dietary changes

**Research and Development:**
- Investigating amyloid-beta and tau-targeted therapies
- Exploring the role of neuroinflammation and the immune system
- Developing biomarkers for early detection and monitoring disease progression
- Testing gene therapy and stem cell approaches

**Further Reading:** [Alzheimer's Association](https://www.alz.org/alzheimers-dementia/what-is-alzheimers)
""",
    "Parkinson's Disease": """
Parkinson's disease is a neurodegenerative disorder that affects predominately dopamine-producing neurons in a specific area of the brain called substantia nigra. Symptoms generally develop slowly over years and include tremors, stiffness, and difficulty with balance and coordination.

**Symptoms:**
- Tremor, mainly at rest (pill-rolling tremor in hands)
- Bradykinesia (slowness of movement)
- Muscle stiffness (rigidity)
- Impaired posture and balance
- Loss of automatic movements (blinking, smiling, swinging arms while walking)
- Speech changes (soft or slurred speech)
- Writing changes (small and crowded handwriting)

**Pathophysiology:**
- Loss of dopaminergic neurons in the substantia nigra
- Accumulation of Lewy bodies (abnormal aggregates of protein) in neurons
- Reduced levels of dopamine in the basal ganglia affecting motor control

**Current Treatments:**
- Dopaminergic medications (levodopa/carbidopa) to replenish dopamine levels
- Dopamine agonists (pramipexole, ropinirole) to mimic dopamine effects
- MAO-B inhibitors (selegiline, rasagiline) to prevent dopamine breakdown
- Deep brain stimulation (DBS) for advanced stages
- Physical, occupational, and speech therapy to manage symptoms

**Research and Development:**
- Developing neuroprotective therapies to slow disease progression
- Investigating alpha-synuclein-targeted treatments
- Exploring the role of gut microbiome and environmental factors
- Advancing gene therapy and stem cell research

**Further Reading:** [Parkinson's Foundation](https://www.parkinson.org/understanding-parkinsons/what-is-parkinsons)
""",
    "Epilepsy": """
Epilepsy is a central nervous system (neurological) disorder in which brain activity becomes abnormal, causing seizures or periods of unusual behavior, sensations, and sometimes loss of awareness.

**Symptoms:**
- Temporary confusion
- A staring spell
- Uncontrollable jerking movements of the arms and legs (convulsions)
- Loss of consciousness or awareness
- Psychic symptoms such as fear, anxiety, or déjà vu

**Types of Seizures:**
- Focal seizures (affecting one part of the brain)
  - Simple focal seizures: No loss of consciousness
  - Complex focal seizures: Altered awareness or consciousness
- Generalized seizures (affecting both sides of the brain)
  - Absence seizures: Brief loss of awareness (staring spells)
  - Tonic-clonic seizures: Convulsions, muscle stiffness, and loss of consciousness
  - Myoclonic seizures: Sudden, brief jerks or twitches
  - Atonic seizures: Sudden loss of muscle tone

**Pathophysiology:**
- Abnormal electrical activity in the brain due to genetic factors, brain injury, or developmental issues
- Imbalance between excitatory and inhibitory neurotransmission

**Current Treatments:**
- Anti-epileptic drugs (AEDs) to control seizures (e.g., phenytoin, valproate, lamotrigine)
- Ketogenic diet (high-fat, low-carbohydrate) for drug-resistant epilepsy
- Vagus nerve stimulation (VNS)
- Responsive neurostimulation (RNS)
- Epilepsy surgery to remove or isolate the seizure focus

**Research and Development:**
- Developing novel AEDs with fewer side effects
- Investigating gene therapy and molecular targets for precision medicine
- Exploring the role of inflammation and immune response in epilepsy
- Advancing brain-computer interface (BCI) technologies for seizure prediction and control

**Further Reading:** [Epilepsy Foundation](https://www.epilepsy.com/learn/about-epilepsy-basics/what-epilepsy)
""",
    "Multiple Sclerosis": """
Multiple sclerosis (MS) is a potentially disabling disease of the brain and spinal cord (central nervous system). In MS, the immune system attacks the protective sheath (myelin) that covers nerve fibers and causes communication problems between your brain and the rest of your body.

**Symptoms:**
- Numbness or weakness in one or more limbs
- Electric-shock sensations with certain neck movements
- Tremor, lack of coordination, or unsteady gait
- Vision problems (blurred vision, double vision, partial or complete vision loss)
- Prolonged double vision
- Slurred speech
- Fatigue
- Dizziness
- Tingling or pain in parts of your body

**Types of MS:**
- Relapsing-remitting MS (RRMS): Periods of new symptoms or relapses followed by periods of remission
- Secondary-progressive MS (SPMS): Initial relapsing-remitting course, followed by progressive worsening
- Primary-progressive MS (PPMS): Gradual onset and steady progression of symptoms without relapses
- Progressive-relapsing MS (PRMS): Steady progression from the beginning with occasional exacerbations

**Pathophysiology:**
- Autoimmune attack on myelin and axons in the central nervous system
- Formation of scar tissue (sclerosis) leading to disrupted nerve signals
- Inflammation and neurodegeneration

**Current Treatments:**
- Disease-modifying therapies (DMTs) to slow progression (e.g., interferons, glatiramer acetate, monoclonal antibodies)
- Corticosteroids to reduce inflammation during relapses
- Physical therapy to improve mobility and strength
- Medications to manage symptoms such as muscle spasms, pain, and bladder dysfunction

**Research and Development:**
- Investigating remyelination and neuroprotection strategies
- Developing advanced imaging techniques for early diagnosis and monitoring
- Exploring the role of gut microbiota in MS
- Testing stem cell therapies and regenerative medicine

**Further Reading:** [National MS Society](https://www.nationalmssociety.org/What-is-MS)
""",
    "Amyotrophic Lateral Sclerosis (ALS)": """
Amyotrophic lateral sclerosis (ALS) is a progressive neurodegenerative disease that affects nerve cells in the brain and the spinal cord. It causes the degeneration of motor neurons, leading to muscle weakness, atrophy, and eventually paralysis.

**Symptoms:**
- Difficulty walking or doing normal daily activities
- Tripping and falling
- Weakness in your leg, feet, or ankles
- Hand weakness or clumsiness
- Slurred speech or trouble swallowing
- Muscle cramps and twitching in your arms, shoulders, and tongue
- Difficulty holding your head up or keeping a good posture
- Uncontrollable laughing or crying (pseudobulbar affect)
- Cognitive and behavioral changes

**Pathophysiology:**
- Progressive degeneration and death of motor neurons in the brain and spinal cord
- Accumulation of protein aggregates (e.g., TDP-43) in neurons
- Involvement of oxidative stress, mitochondrial dysfunction, and excitotoxicity

**Current Treatments:**
- Riluzole and edaravone to slow disease progression
- Medications to manage symptoms such as muscle cramps, spasticity, and pain
- Non-invasive ventilation to assist breathing
- Nutritional support to ensure adequate calorie intake
- Physical and occupational therapy to maintain mobility and independence
- Speech therapy to address communication difficulties

**Research and Development:**
- Investigating genetic and molecular mechanisms of ALS
- Developing gene therapy and RNA-based treatments
- Exploring the role of neuroinflammation and immune modulation
- Testing novel neuroprotective and regenerative approaches

**Further Reading:** [ALS Association](https://www.als.org/understanding-als/what-is-als)
""",
    "Huntington's Disease": """
Huntington's disease is a genetic disorder that causes the progressive breakdown of nerve cells in the brain. It deteriorates a person's physical and mental abilities, usually during their prime working years, and has no cure.

**Symptoms:**
- Involuntary jerking or writhing movements (chorea)
- Muscle problems, such as rigidity or muscle contracture (dystonia)
- Slow or abnormal eye movements
- Impaired gait, posture, and balance
- Difficulty with speech or swallowing
- Cognitive impairments such as difficulty organizing, prioritizing, or focusing on tasks
- Lack of flexibility or the tendency to get stuck on a thought, behavior, or action (perseveration)
- Lack of impulse control that can result in outbursts, acting without thinking, and sexual promiscuity
- Lack of awareness of one’s own behaviors and abilities
- Slowness in processing thoughts or ''finding'' words
- Difficulty in learning new information
- Depression, anxiety, and other psychiatric disorders

**Pathophysiology:**
- Caused by a genetic mutation in the HTT gene, leading to an abnormal expansion of the CAG trinucleotide repeat
- Accumulation of mutant huntingtin protein in neurons
- Progressive degeneration of neurons in the basal ganglia and cortex

**Current Treatments:**
- Medications to manage symptoms such as chorea (tetrabenazine, deutetrabenazine)
- Antipsychotic medications to manage psychiatric symptoms
- Antidepressants to treat depression and anxiety
- Physical therapy to improve strength, flexibility, and balance
- Speech therapy to address communication and swallowing difficulties
- Nutritional support to ensure adequate calorie intake

**Research and Development:**
- Investigating gene-silencing approaches (e.g., antisense oligonucleotides, RNA interference)
- Exploring neuroprotective strategies and regenerative medicine
- Testing novel therapies targeting mutant huntingtin protein
- Developing biomarkers for early diagnosis and monitoring disease progression

**Further Reading:** [Huntington's Disease Society of America](https://hdsa.org/what-is-hd/)
""",
    "Migraine": """
Migraine is a neurological condition that can cause multiple symptoms. It's frequently characterized by intense, debilitating headaches. Symptoms may include nausea, vomiting, difficulty speaking, numbness or tingling, and sensitivity to light and sound.

**Symptoms:**
- Severe, throbbing or pulsating headache
- Pain usually on one side of the head, but often on both sides
- Sensitivity to light, sound, and sometimes smell and touch
- Nausea and vomiting
- Blurred vision
- Lightheadedness, sometimes followed by fainting
- Aura (visual disturbances such as flashing lights or blind spots)
- Tingling or numbness in the face or extremities

**Pathophysiology:**
- Activation of the trigeminovascular system leading to neurogenic inflammation
- Involvement of neurotransmitters such as serotonin and calcitonin gene-related peptide (CGRP)
- Genetic and environmental factors contributing to susceptibility

**Current Treatments:**
- Acute treatments to relieve symptoms (e.g., triptans, NSAIDs, antiemetics)
- Preventive treatments to reduce the frequency and severity of migraines (e.g., beta-blockers, anticonvulsants, CGRP inhibitors)
- Lifestyle changes to avoid triggers (e.g., stress management, regular sleep patterns, dietary modifications)
- Alternative therapies such as acupuncture, biofeedback, and cognitive-behavioral therapy

**Research and Development:**
- Developing novel acute and preventive treatments targeting CGRP and other pathways
- Investigating the role of the gut-brain axis in migraine
- Exploring the genetic basis of migraine susceptibility
- Advancing neuroimaging techniques to study migraine pathophysiology

**Further Reading:** [Migraine Research Foundation](https://migraineresearchfoundation.org/about-migraine/migraine-facts/)
""",
    "Autism Spectrum Disorder (ASD)": """
Autism Spectrum Disorder (ASD) is a developmental disorder that affects communication and behavior. It encompasses a range of conditions characterized by challenges with social skills, repetitive behaviors, speech, and nonverbal communication.

**Symptoms:**
- Difficulty with communication and interaction with other people
- Restricted interests and repetitive behaviors
- Symptoms that affect the individual’s ability to function in school, work, and other areas of life
- Avoidance of eye contact and preference for being alone
- Difficulty understanding other people’s feelings or talking about their own
- Delay in speech and language skills
- Repeating words or phrases over and over (echolalia)
- Unusual reactions to the way things sound, smell, taste, look, or feel

**Pathophysiology:**
- Genetic mutations and environmental factors affecting brain development
- Abnormalities in brain structure and function
- Imbalance in excitatory and inhibitory neurotransmission

**Current Treatments:**
- Behavioral therapy to improve social skills and communication
- Speech therapy to address language difficulties
- Occupational therapy to help with daily living skills
- Medications to manage symptoms like anxiety, depression, and hyperactivity
- Educational interventions tailored to individual needs

**Research and Development:**
- Investigating the genetic and molecular basis of ASD
- Exploring the role of the gut-brain axis and immune system in ASD
- Developing novel therapies targeting synaptic function and neuroinflammation
- Advancing early diagnosis and intervention strategies

**Further Reading:** [Autism Speaks](https://www.autismspeaks.org/what-autism)
""",
    "Attention-Deficit/Hyperactivity Disorder (ADHD)": """
Attention-Deficit/Hyperactivity Disorder (ADHD) is a neurodevelopmental disorder characterized by inattention, hyperactivity, and impulsiveness. It is commonly diagnosed in children but can persist into adulthood.

**Symptoms:**
- Inattention (difficulty staying focused, disorganized, easily distracted)
- Hyperactivity (excessive movement, fidgeting, difficulty sitting still)
- Impulsiveness (hasty actions without thinking, interrupting others, difficulty waiting)

**Types of ADHD:**
- Predominantly inattentive presentation
- Predominantly hyperactive-impulsive presentation
- Combined presentation

**Pathophysiology:**
- Imbalance in neurotransmitters such as dopamine and norepinephrine
- Differences in brain structure and function, particularly in areas related to attention and executive function
- Genetic and environmental factors contributing to susceptibility

**Current Treatments:**
- Stimulant medications (e.g., methylphenidate, amphetamines) to improve attention and focus
- Non-stimulant medications (e.g., atomoxetine, guanfacine) for those who do not respond to stimulants
- Behavioral therapy to develop coping strategies and improve behavior
- Educational interventions to support learning and academic performance
- Parent training programs to manage behavior and support development

**Research and Development:**
- Investigating the genetic and neurobiological basis of ADHD
- Exploring the role of diet, sleep, and lifestyle factors in ADHD management
- Developing novel pharmacological treatments with fewer side effects
- Advancing early diagnosis and personalized treatment approaches

**Further Reading:** [CHADD](https://chadd.org/about-adhd/overview/)
""",
    "Stroke": """
A stroke occurs when the blood supply to part of your brain is interrupted or reduced, preventing brain tissue from getting oxygen and nutrients. Brain cells begin to die in minutes.

**Types of Stroke:**
- Ischemic stroke: Caused by a blockage in an artery supplying blood to the brain.
- Hemorrhagic stroke: Caused by a burst blood vessel in the brain.
- Transient ischemic attack (TIA): A temporary period of symptoms similar to those of a stroke.

**Symptoms:**
- Sudden numbness or weakness in the face, arm, or leg, especially on one side of the body
- Sudden confusion, trouble speaking, or difficulty understanding speech
- Sudden trouble seeing in one or both eyes
- Sudden trouble walking, dizziness, loss of balance, or lack of coordination
- Sudden severe headache with no known cause

**Pathophysiology:**
- Ischemic stroke: Blood clot blocks or narrows an artery leading to the brain, causing reduced blood flow and oxygen deprivation.
- Hemorrhagic stroke: Blood vessel ruptures, causing bleeding in or around the brain, leading to increased pressure and damage to brain cells.

**Current Treatments:**
- Emergency treatments to restore blood flow (e.g., clot-busting medications, mechanical thrombectomy)
- Medications to prevent further strokes (e.g., antiplatelets, anticoagulants)
- Rehabilitation therapies (e.g., physical therapy, occupational therapy, speech therapy) to regain lost functions
- Lifestyle changes to manage risk factors (e.g., controlling blood pressure, quitting smoking, maintaining a healthy diet)

**Research and Development:**
- Investigating neuroprotective agents to minimize brain damage during stroke
- Exploring stem cell therapy and regenerative medicine to repair brain damage
- Developing advanced imaging techniques for early diagnosis and treatment planning
- Testing novel therapies to improve recovery and reduce disability

**Further Reading:** [American Stroke Association](https://www.stroke.org/en/about-stroke)
""",
    "Amyotrophic Lateral Sclerosis (ALS)": """
Amyotrophic lateral sclerosis (ALS) is a progressive neurodegenerative disease that affects nerve cells in the brain and the spinal cord. It causes the degeneration of motor neurons, leading to muscle weakness, atrophy, and eventually paralysis.

**Symptoms:**
- Difficulty walking or doing normal daily activities
- Tripping and falling
- Weakness in your leg, feet, or ankles
- Hand weakness or clumsiness
- Slurred speech or trouble swallowing
- Muscle cramps and twitching in your arms, shoulders, and tongue
- Difficulty holding your head up or keeping a good posture
- Uncontrollable laughing or crying (pseudobulbar affect)
- Cognitive and behavioral changes

**Pathophysiology:**
- Progressive degeneration and death of motor neurons in the brain and spinal cord
- Accumulation of protein aggregates (e.g., TDP-43) in neurons
- Involvement of oxidative stress, mitochondrial dysfunction, and excitotoxicity

**Current Treatments:**
- Riluzole and edaravone to slow disease progression
- Medications to manage symptoms such as muscle cramps, spasticity, and pain
- Non-invasive ventilation to assist breathing
- Nutritional support to ensure adequate calorie intake
- Physical and occupational therapy to maintain mobility and independence
- Speech therapy to address communication difficulties

**Research and Development:**
- Investigating genetic and molecular mechanisms of ALS
- Developing gene therapy and RNA-based treatments
- Exploring the role of neuroinflammation and immune modulation
- Testing novel neuroprotective and regenerative approaches

**Further Reading:** [ALS Association](https://www.als.org/understanding-als/what-is-als)
""",
    "Cerebral Palsy": """
Cerebral palsy is a group of disorders that affect a person's ability to move and maintain balance and posture. It is the most common motor disability in childhood.

**Symptoms:**
- Variations in muscle tone (too stiff or too floppy)
- Spasticity (stiff muscles and exaggerated reflexes)
- Ataxia (lack of muscle coordination)
- Tremors or involuntary movements
- Delays in reaching motor skills milestones
- Difficulty walking
- Excessive drooling or problems with swallowing
- Difficulty with precise movements (e.g., picking up objects)

**Pathophysiology:**
- Caused by abnormal brain development or damage to the developing brain
- Factors include genetic mutations, maternal infections, fetal stroke, lack of oxygen during birth, traumatic head injury, and brain infections

**Current Treatments:**
- Medications to manage symptoms such as muscle spasticity and seizures
- Physical therapy to improve muscle strength and coordination
- Occupational therapy to assist with daily living activities
- Speech therapy to improve communication and swallowing difficulties
- Orthopedic surgery to correct bone and joint deformities
- Assistive devices (e.g., braces, wheelchairs) to improve mobility

**Research and Development:**
- Investigating early interventions to improve outcomes
- Exploring stem cell therapy and regenerative medicine
- Developing advanced neuroimaging techniques for diagnosis and monitoring
- Testing novel pharmacological treatments to manage symptoms

**Further Reading:** [Cerebral Palsy Foundation](https://www.yourcpf.org/what-is-cerebral-palsy/)
""",
    "Schizophrenia": """
Schizophrenia is a chronic brain disorder that affects less than one percent of the U.S. population. When schizophrenia is active, symptoms can include delusions, hallucinations, trouble with thinking and concentration, and lack of motivation.

**Symptoms:**
- Delusions (false beliefs)
- Hallucinations (seeing or hearing things that are not there)
- Disorganized thinking (thoughts may be jumbled or make no sense)
- Extremely disorganized or abnormal motor behavior
- Negative symptoms (reduced ability to function normally)

**Pathophysiology:**
- Involves imbalances in neurotransmitters such as dopamine and glutamate
- Genetic and environmental factors contributing to the risk
- Structural and functional abnormalities in the brain, particularly in the prefrontal cortex and hippocampus

**Current Treatments:**
- Antipsychotic medications to manage symptoms (e.g., risperidone, olanzapine)
- Psychosocial interventions (e.g., cognitive-behavioral therapy, social skills training)
- Supportive therapy to help manage daily life and improve quality of life
- Rehabilitation programs to improve social and vocational skills
- Coordinated specialty care programs for early intervention and support

**Research and Development:**
- Investigating the genetic and neurobiological basis of schizophrenia
- Exploring the role of inflammation and immune system in schizophrenia
- Developing novel antipsychotic medications with fewer side effects
- Testing advanced neuroimaging techniques for early diagnosis and monitoring

**Further Reading:** [National Institute of Mental Health](https://www.nimh.nih.gov/health/topics/schizophrenia)
""",
    "Bipolar Disorder": """
Bipolar disorder, formerly called manic-depressive illness or manic depression, is a mental illness that causes unusual shifts in mood, energy, activity levels, concentration, and the ability to carry out day-to-day tasks.

**Symptoms:**
- Manic episodes (elevated mood, increased activity, inflated self-esteem, decreased need for sleep)
- Depressive episodes (low mood, fatigue, feelings of worthlessness, difficulty concentrating, changes in sleep and appetite)
- Hypomanic episodes (less severe manic symptoms)
- Mixed episodes (symptoms of both mania and depression)

**Types of Bipolar Disorder:**
- Bipolar I Disorder: Characterized by manic episodes that last at least 7 days, or by manic symptoms that are so severe that immediate hospital care is needed.
- Bipolar II Disorder: Defined by a pattern of depressive episodes and hypomanic episodes, but not the full-blown manic episodes that are typical of Bipolar I Disorder.
- Cyclothymic Disorder (also called Cyclothymia): Periods of hypomanic symptoms as well as periods of depressive symptoms lasting for at least 2 years (1 year in children and adolescents).

**Pathophysiology:**
- Involves imbalances in neurotransmitters such as dopamine, serotonin, and norepinephrine
- Genetic predisposition and environmental triggers
- Abnormalities in brain circuits regulating mood and behavior

**Current Treatments:**
- Mood stabilizers (e.g., lithium, valproate) to manage mood swings
- Antipsychotic medications (e.g., quetiapine, olanzapine) for manic episodes
- Antidepressants for depressive episodes
- Psychotherapy (e.g., cognitive-behavioral therapy, family-focused therapy)
- Lifestyle changes (e.g., regular sleep patterns, stress management)

**Research and Development:**
- Investigating the genetic and molecular basis of bipolar disorder
- Exploring the role of inflammation and immune system in mood disorders
- Developing novel pharmacological treatments with fewer side effects
- Advancing early diagnosis and personalized treatment approaches

**Further Reading:** [National Institute of Mental Health](https://www.nimh.nih.gov/health/topics/bipolar-disorder)
""",
    "Obsessive-Compulsive Disorder (OCD)": """
Obsessive-Compulsive Disorder (OCD) is a chronic mental health condition characterized by uncontrollable, reoccurring thoughts (obsessions) and behaviors (compulsions) that the person feels the urge to repeat over and over.

**Symptoms:**
- Obsessions (repeated, persistent, and unwanted thoughts, urges, or images)
- Compulsions (repetitive behaviors or mental acts aimed at reducing anxiety related to obsessions)
- Common themes include fear of contamination, need for symmetry, aggressive thoughts, and fear of harm

**Pathophysiology:**
- Involves imbalances in neurotransmitters such as serotonin
- Abnormalities in brain circuits regulating behavior and emotion
- Genetic and environmental factors contributing to risk

**Current Treatments:**
- Selective serotonin reuptake inhibitors (SSRIs) to reduce symptoms (e.g., fluoxetine, sertraline)
- Cognitive-behavioral therapy (CBT) with exposure and response prevention (ERP)
- Deep brain stimulation (DBS) for severe, treatment-resistant cases
- Psychosocial interventions to improve coping and functioning

**Research and Development:**
- Investigating the genetic and neurobiological basis of OCD
- Exploring novel pharmacological treatments targeting different neurotransmitter systems
- Developing advanced neuroimaging techniques for diagnosis and monitoring
- Testing novel therapeutic approaches such as transcranial magnetic stimulation (TMS)

**Further Reading:** [International OCD Foundation](https://iocdf.org/about-ocd/)
""",
    "Tourette Syndrome": """
Tourette Syndrome (TS) is a neurological disorder characterized by repetitive, stereotyped, involuntary movements and vocalizations called tics. It often begins in childhood.

**Symptoms:**
- Motor tics (sudden, brief, repetitive movements such as blinking, shrugging, or facial grimacing)
- Vocal tics (repetitive sounds such as throat clearing, grunting, or shouting)
- Tics can be simple (involving a few parts of the body) or complex (involving several parts of the body or words/phrases)

**Pathophysiology:**
- Involves abnormalities in brain circuits regulating movement and behavior
- Genetic and environmental factors contributing to risk
- Imbalances in neurotransmitters such as dopamine

**Current Treatments:**
- Medications to manage tics (e.g., antipsychotics, alpha-adrenergic agonists)
- Behavioral therapy (e.g., habit reversal training, comprehensive behavioral intervention for tics)
- Supportive therapies to address associated conditions such as ADHD, OCD, and anxiety
- Educational interventions to support learning and social integration

**Research and Development:**
- Investigating the genetic and neurobiological basis of TS
- Exploring novel pharmacological treatments targeting different neurotransmitter systems
- Developing advanced neuroimaging techniques for diagnosis and monitoring
- Testing novel therapeutic approaches such as deep brain stimulation (DBS)

**Further Reading:** [Tourette Association of America](https://tourette.org/about-tourette/overview/)
""",
    "Restless Legs Syndrome (RLS)": """
Restless Legs Syndrome (RLS), also known as Willis-Ekbom Disease, is a neurological disorder characterized by an irresistible urge to move the legs, often accompanied by uncomfortable sensations.

**Symptoms:**
- Unpleasant sensations in the legs (itching, tingling, pulling, or crawling)
- Strong urge to move the legs to relieve discomfort
- Symptoms typically occur in the evening or at night
- Difficulty falling asleep or staying asleep due to symptoms
- Daytime fatigue and impaired functioning

**Pathophysiology:**
- Involves abnormalities in dopamine regulation in the brain
- Genetic predisposition and environmental triggers
- Possible involvement of iron deficiency and peripheral neuropathy

**Current Treatments:**
- Medications to manage symptoms (e.g., dopamine agonists, anticonvulsants, opioids)
- Iron supplements for individuals with low iron levels
- Lifestyle changes to improve sleep hygiene and reduce symptoms (e.g., regular exercise, avoiding caffeine)
- Behavioral therapy to address sleep disturbances and anxiety

**Research and Development:**
- Investigating the genetic and neurobiological basis of RLS
- Exploring novel pharmacological treatments targeting dopamine and other neurotransmitter systems
- Developing advanced diagnostic criteria and tools for early detection
- Testing novel therapeutic approaches such as transcranial magnetic stimulation (TMS)

**Further Reading:** [Restless Legs Syndrome Foundation](https://www.rls.org/understanding-rls)
""",
    "Myasthenia Gravis": """
Myasthenia Gravis (MG) is a chronic autoimmune neuromuscular disorder characterized by weakness and rapid fatigue of voluntary muscles. It occurs when the immune system mistakenly attacks the connections between nerves and muscles.

**Symptoms:**
- Weakness in the arm and leg muscles
- Double vision and drooping eyelids (ptosis)
- Difficulties with speech, chewing, swallowing, and breathing
- Symptoms may fluctuate in intensity and improve with rest

**Pathophysiology:**
- Autoimmune attack on acetylcholine receptors at the neuromuscular junction
- Reduction in the number of functional receptors leading to impaired muscle contraction
- Genetic and environmental factors contributing to risk

**Current Treatments:**
- Medications to improve neuromuscular transmission (e.g., acetylcholinesterase inhibitors)
- Immunosuppressive therapies to reduce the autoimmune response (e.g., corticosteroids, azathioprine)
- Plasmapheresis and intravenous immunoglobulin (IVIG) for severe cases
- Thymectomy (surgical removal of the thymus gland) in some patients

**Research and Development:**
- Investigating the genetic and immunological basis of MG
- Exploring novel immunomodulatory treatments and targeted therapies
- Developing advanced diagnostic criteria and tools for early detection
- Testing novel therapeutic approaches such as monoclonal antibodies

**Further Reading:** [Myasthenia Gravis Foundation of America](https://myasthenia.org/What-is-MG)
""",
    "Narcolepsy": """
Narcolepsy is a chronic neurological disorder that affects the brain's ability to regulate sleep-wake cycles. It is characterized by excessive daytime sleepiness, cataplexy (sudden muscle weakness), sleep paralysis, and hallucinations.

**Symptoms:**
- Excessive daytime sleepiness (EDS) leading to sudden sleep attacks
- Cataplexy triggered by strong emotions such as laughter or surprise
- Sleep paralysis (temporary inability to move or speak while falling asleep or waking up)
- Hypnagogic hallucinations (vivid, dream-like experiences while falling asleep)
- Disrupted nighttime sleep

**Pathophysiology:**
- Loss of hypocretin (orexin) neurons in the hypothalamus, leading to impaired regulation of sleep-wake cycles
- Genetic predisposition and environmental triggers
- Abnormalities in REM sleep regulation

**Current Treatments:**
- Stimulant medications to improve wakefulness (e.g., modafinil, amphetamines)
- Medications to manage cataplexy and REM sleep abnormalities (e.g., sodium oxybate, SSRIs, SNRIs)
- Lifestyle modifications to improve sleep hygiene and manage symptoms
- Scheduled naps to reduce daytime sleepiness

**Research and Development:**
- Investigating the genetic and neurobiological basis of narcolepsy
- Exploring novel pharmacological treatments targeting hypocretin and other neurotransmitter systems
- Developing advanced diagnostic criteria and tools for early detection
- Testing novel therapeutic approaches such as immunotherapy

**Further Reading:** [Narcolepsy Network](https://narcolepsynetwork.org/about-narcolepsy/)
""",
    "Guillain-Barre Syndrome (GBS)": """
Guillain-Barre Syndrome (GBS) is a rare neurological disorder in which the body's immune system mistakenly attacks the peripheral nerves, leading to muscle weakness, numbness, and paralysis.

**Symptoms:**
- Rapid onset of muscle weakness, starting in the legs and spreading to the upper body
- Tingling or numbness in the extremities
- Difficulty with walking, coordination, and balance
- Severe cases can lead to paralysis and respiratory failure

**Pathophysiology:**
- Autoimmune attack on the myelin sheath or axons of peripheral nerves
- Often triggered by infections such as respiratory or gastrointestinal infections
- Involvement of both humoral and cellular immune responses

**Current Treatments:**
- Intravenous immunoglobulin (IVIG) to modulate the immune response
- Plasmapheresis (plasma exchange) to remove antibodies from the blood
- Physical therapy to support recovery and improve mobility
- Supportive care to manage complications such as respiratory support and pain management

**Research and Development:**
- Investigating the genetic and immunological basis of GBS
- Exploring novel immunomodulatory treatments and targeted therapies
- Developing advanced diagnostic criteria and tools for early detection
- Testing novel therapeutic approaches such as complement inhibitors

**Further Reading:** [Guillain-Barre Syndrome Foundation](https://www.gbs-cidp.org/about-gbs/)
"""
}

# Title and introductory text
st.title("Neurological Diseases Information")
st.write("Select a neurological disease from the dropdown menu below to learn more about it.")

# Dropdown menu for selecting a neurological disease
disease = st.selectbox(
    "Choose a neurological disease:",
    list(disease_info.keys())
)

# Display information about the selected disease
st.subheader(disease)
st.markdown(disease_info[disease])


st.write("---")
st.markdown("""
    <div style='display: flex; justify-content: space-around; padding: 20px; background-color: white;'>
        <div style='background-color: white; padding: 10px;'>
            <p><a href='#'>FAQ</a></p>
        </div>
        <div style='background-color: white; padding: 10px;'>
            <p><a href='#'>User Guide</a></p>
        </div>
        <div style='background-color: white; padding: 10px;'>
            <p><a href='#'>Privacy Policy</a></p>
        </div>
        <div style='background-color: white; padding: 10px;'>
            <p><a href='#'>Terms of Service</a></p>
        </div>
    </div>
    <hr style='margin: 20px 0; padding: 0;'>
    <div style='text-align: center; padding: 10px;'>
        &copy; 2024 Drug Development App. All rights reserved. | Contact: info@drugdevapp.com
    </div>
    """, unsafe_allow_html=True)
