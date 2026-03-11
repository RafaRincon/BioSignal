# BioSignal Frontend

React + Vite + Express. Vive en `frontend/` dentro del repo principal.

## Estructura

```
frontend/
├── index.html
├── vite.config.js
├── package.json
├── server/
│   ├── index.js            ← Express (puerto 3001)
│   └── routes/
│       └── reports.js      ← GET /api/reports/:disease
└── src/
    ├── main.jsx
    ├── App.jsx
    ├── styles/
    │   └── globals.css     ← tokens del sistema de diseño
    ├── services/
    │   └── api.js          ← ÚNICO punto de contacto con backend
    ├── hooks/
    │   └── useReport.js    ← fetch + loading + error
    ├── utils/
    │   └── format.js       ← funciones puras de formato
    └── components/
        ├── Sidebar.jsx
        └── overview/
            ├── OverviewTab.jsx         ← orquestador
            ├── PageHeader.jsx
            ├── PipelineStatus.jsx
            ├── SummaryHero.jsx
            ├── ConsensusGenes.jsx
            ├── TopPathways.jsx
            ├── TherapeuticTargets.jsx
            ├── KeyMechanisms.jsx
            ├── NovelHypotheses.jsx
            └── BiomarkerCandidates.jsx
```

## Instalación

```bash
# Desde la raíz del repo
cd frontend
npm install
npm install cors   # dependencia del servidor Express
```

## Desarrollo

Terminal 1 — backend Express:
```bash
node server/index.js
```

Terminal 2 — frontend Vite:
```bash
npm run dev
```

Abre http://localhost:5173

## Requisito de datos

El backend lee reportes desde:
```
data/reports/{disease}/report_{disease}_{date}.json
```

El JSON de Alzheimer ya está en:
```
data/reports/alzheimer/report_alzheimer_20260310.json
```

## Producción

```bash
npm run build
NODE_ENV=production node server/index.js
```

## Agregar un nuevo tab

1. Crea `src/components/{tab}/` con sus sub-componentes
2. Importa en `src/App.jsx` y agrega la condición de render
3. Agrega el item al array `NAV_*` en `src/components/Sidebar.jsx`
4. Si necesita un nuevo endpoint, agrégalo en `src/services/api.js` y crea la ruta en `server/routes/`
